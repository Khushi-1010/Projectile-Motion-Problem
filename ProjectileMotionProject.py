"""
====================================================================================================
SCRIPT NAME: ProjectileMotionProject.py
AUTHOR: Khushi Piparava
INITIATED: November 18, 2024
LAST REVISION: November 27, 2024
====================================================================================================
SCRIPT DESCRIPTION:
This script simulates the trajectory of the given preojectile in ENZ coordinates.
====================================================================================================
USER-DEFINED FUNCTIONS:
constants()|This function loads all of the constants used to solve the free fall problem.
----------------------------------------------------------------------------------------------------
drag()|This function calculates calculates the drag force in the ENZ frame based on the vehicle's 
position, velocity, and aerodynamic properties.
----------------------------------------------------------------------------------------------------
gravity()|This function calculates the acceleration due to gravity in ENZ frame.
----------------------------------------------------------------------------------------------------
projectile_eom()|This function calculates the equations of motion for the free fall problem in state
form.
----------------------------------------------------------------------------------------------------
plot_position()|This function plots the vehicle's displacements as a function of time.
----------------------------------------------------------------------------------------------------
plot_velocity()|This function plots the vehicle's speeds as a function of time.
----------------------------------------------------------------------------------------------------
plot_acceleration()|This function plots the vehicle's inertial acceleration as a function of time.
----------------------------------------------------------------------------------------------------
plot_altitude()|This function plots the vehicle's linear range as a function of time.
----------------------------------------------------------------------------------------------------
plot_angles()|This function plots the vehicle's Azimuth and Elevation angles as a function of time.
----------------------------------------------------------------------------------------------------
plot_pressure()|This function plots the Dynamic Pressure as a function of time.
----------------------------------------------------------------------------------------------------
run_simulation()|This function runs the simulation
====================================================================================================
ABBREVIATIONS:
'GS' = "Ground Station".
----------------------------------------------------------------------------------------------------
'CG' = "Center of Gravity of the projectile".
----------------------------------------------------------------------------------------------------
'ENZ' = "East-North-Zenith"
----------------------------------------------------------------------------------------------------
'WRT' = "with respect to".
====================================================================================================
ADDITIONAL COMMENTS:
None.
====================================================================================================
PERMISSION:
Any use of this code either in part or in full must first be approved by Khushi Piparava, Student in 
the Department of Mechanical and Aerospace Engineering at The University of Texas at Arlington.  
For permission to use this code, Khushi may be contacted at krp3505@mavs.uta.edu.
====================================================================================================
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import time

def constants():
    # Dictionary of constants
    C = {
        'u': 398600.435436E9, # [m^3/s^2]Gravitational parameter of Earth
        'Re': 6378.137E3,   # [m]Mean equatorial radius of Earth
        'we': 2 * np.pi / 86164.1,  # [rad/s]Rotational speed of Earth
        'g': 9.80665,   # [m/s^2]Standard Acceleration due to Gravity
        'GS': {
            'Lat': (32 + 44 / 60 + 52 / 3600) * np.pi / 180, # GS Latitude
            'Long': -(97 + 5 / 60 + 34 / 3600) * np.pi / 180, # GS Longitude
            'hgs': 185, # GS Altitude above mean equator
        },
        'CG': {
            'Cd': 0.82, # []Drag coefficient
            'vcggs': 827, # [m/s]Speed of projectile
            'Azimuth': 107 * np.pi / 180, # [rad] Azimuth of projectile
            'Elevation': 60 * np.pi / 180, # [rad] Elevation of projectile
            'mcg': 35, # [kg]Vehicle Mass
            'Ar': 6.00225 * 10**-3 * np.pi, # [m^2]Vehicle reference area
            'Rcggs': np.array([0, 0, 0])    # [m]Initial CG position relative to the GS in ENZ
        }
    }

    C['GS']['We'] = C['we'] * np.array([0, np.cos(C['GS']['Lat']), np.sin(C['GS']['Lat'])])
    # [rad/s]Rotational velocity of the Earth in ENZ
    
    C['GS']['Rgse'] = (C['Re'] + C['GS']['hgs']) * np.array([0, 0, 1])
    # [m]GS position wrt Earth in ENZ

    C['GS']['Vgse'] = np.cross(C['GS']['We'], C['GS']['Rgse'])
    # [m/s]Ground station velocity WRT the Earth in ENZ

    C['GS']['Agse'] = np.cross(C['GS']['We'], C['GS']['Vgse'])
    # [m/s^2]Ground station acceleration WRT the Earth ENZ

    C['CG']['Vcggs'] = C['CG']['vcggs'] * np.array([
        np.cos(C['CG']['Elevation']) * np.sin(C['CG']['Azimuth']),
        np.cos(C['CG']['Elevation']) * np.cos(C['CG']['Azimuth']),
        np.sin(C['CG']['Elevation'])
    ]) # [m/s]Initial CG velocity relative to the GS in ENZ

    C['CG']['Scggs'] = np.hstack((C['CG']['Rcggs'], C['CG']['Vcggs']))
    # [m,m/s]Initial CG state relative to the GS in ENZ
    return C

def drag(rcge, Vinf, C):
    vinf = np.linalg.norm(Vinf) # [m/s]True air speed

    hcg = rcge - C['Re'] # [m]Altitude above mean equator
    
    _, _, rho = standard_atmosphere(hcg) # [kg/m^3]Atmospheric density
    
    Fd = -0.5 * C['CG']['Cd'] * rho * C['CG']['Ar'] * vinf * Vinf # [N]Drag Force in ENZ
    
    return Fd

def gravity(Rcge, rcge, C):
    g = -C['u'] * Rcge / (rcge**3) 
    # [m/s^2] Acceleration due to gravity in ENZ
    
    return g

def standard_atmosphere(hcg):
    g = 9.80665
    R = 287
    Re = 6378137
    h = Re * hcg / (Re + hcg)
    
    if h <= 11000:
        a = (216.66 - 288.16) / 11000
        T = 288.16 + a * h
        P = 101325 * (T / 288.16)**(-g / (a * R))
        D = P / (R * T)
    elif h <= 25000:
        T = 216.66
        P = 22650.1684742737 * np.exp(-g / (R * T) * (h - 11000))
        D = P / (R * T)
    elif h <= 47000:
        a = (282.66 - 216.66) / (47000 - 25000)
        T = 216.66 + a * (h - 25000)
        P = 2493.58245271879 * (T / 216.66)**(-g / (a * R))
        D = P / (R * T)
    elif h <= 53000:
        T = 282.66
        P = 120.879682128688 * np.exp(-g / (R * T) * (h - 47000))
        D = P / (R * T)
    elif h <= 79000:
        a = (165.66 - 282.66) / (79000 - 53000)
        T = 282.66 + a * (h - 53000)
        P = 58.5554504138705 * (T / 282.66)**(-g / (a * R))
        D = P / (R * T)
    elif h <= 90000:
        T = 165.66
        P = 1.01573256565262 * np.exp(-g / (R * T) * (h - 79000))
        D = P / (R * T)
    elif h <= 105000:
        a = (225.66 - 165.66) / (105000 - 90000)
        T = 165.66 + a * (h - 90000)
        P = 0.105215646463089 * (T / 165.66)**(-g / (a * R))
        D = P / (R * T)
    elif h <= 500000:
        a = (4 - 225.66) / (500000 - 105000)
        T = 225.66 + a * (h - 105000)
        P = 0.00751891790761519 * (T / 225.66)**(-g / (a * R))
        D = P / (R * T)
    else:
        T = 4
        P = 0
        D = 0
    return T, P, D

def projectile_eom(t, S, C):
    Rcggs = S[:3]
    # [m]Vehicle position relative to the ground station in ENZ

    Vcggs = S[3:6]
    #[m/s]Vehicle velocity relative to the ground station in ENZ

    Rcge = C['GS']['Rgse'] + Rcggs
    # [m]Absolute vehicle position WRT the Earth in ENZ

    rcge = np.linalg.norm(Rcge)
    # [m]Vehicle range WRT the Earth
    
    Vcge = C['GS']['Vgse'] + np.cross(C['GS']['We'], Rcggs) + Vcggs
    # [m/s]Absolute vehicle velocity WRT the Earth in ENZ
    
    Vatm = np.cross(C['GS']['We'], Rcge)
    # [m/s]Atmosphere velocity WRT the Earth in ENZ
    
    Vinf = Vcge - Vatm
    # [m/s]True air velocity
    
    Fd = drag(rcge, Vinf, C)
    # [N]Drag force in ENZ
    
    g = gravity(Rcge, rcge, C)
    # [m/s^2]Acceleration due to gravity

    dSdt = np.zeros(6)
    # []Allocates memory for the state vector derivative
    
    dSdt[:3] = Vcggs
    # [m/s]Vehicle velocity relative to the ground station
    
    dSdt[3:] = Fd / C['CG']['mcg'] + g - C['GS']['Agse'] - np.cross(C['GS']['We'], np.cross(C['GS']['We'], Rcggs)) - 2 * np.cross(C['GS']['We'], Vcggs)
    # [m/s^2]Vehicle acceleration relative to the ground station
    
    return dSdt

def plot_position(t, Rcggs):
    
    titles = ['East', 'North', 'Zenith', 'Linear']
    ylimits = [[0, 7000], [-2500, 0], [0, 6000], [0, 8000]]
    yticks = [
        np.arange(0, 7001, 700),
        np.arange(-2500, 1, 250),
        np.arange(0, 6001, 600),
        np.arange(0, 8001, 800),
    ]
    
    # Validate input dimensions
    if Rcggs.shape[1] != 3:
        raise ValueError("Rcggs must have exactly 3 columns (x, y, z positions).")

    # Create 2x2 subplots
    fig, ax = plt.subplots(2, 2, sharex=True)
    ax = ax.flatten()  # Flatten to access subplots as a 1D array
    #fig.suptitle('Position Analysis', fontsize=16)

    # Plot each coordinate position over time
    for i in range(3):
        ax[i].plot(t, Rcggs[:, i], 'k', linestyle='None', markersize=8, marker='.')
        ax[i].set_title(titles[i], fontsize=16, fontweight='bold', fontname='Arial')
        ax[i].set_ylabel('Displacement (m)',fontsize=12, fontweight='bold', fontname='Arial')
        ax[i].set_ylim(ylimits[i])
        ax[i].set_yticks(yticks[i])
        ax[i].grid()
        
    # Plot total displacement (norm)
    total_displacement = np.linalg.norm(Rcggs, axis=1)
    ax[3].plot(t, total_displacement, 'k', linestyle='None', markersize=8, marker='.')
    ax[3].set_title(titles[3], fontsize=16, fontweight='bold', fontname='Arial')
    ax[3].set_ylabel('Displacement (m)', fontsize=12, fontweight='bold', fontname='Arial')
    ax[3].set_ylim(ylimits[3])
    ax[3].set_yticks(yticks[3])
    ax[3].grid()
    
    # Set x-axis label on the bottom plots
    ax[2].set_xlabel('Time (s)', fontsize=12, fontweight='bold', fontname='Arial')
    ax[3].set_xlabel('Time (s)', fontsize=12, fontweight='bold', fontname='Arial')
    ax[2].set_xlim([0, 70])
    ax[2].set_xticks(np.linspace(0, 70, 7))
    
    fig.tight_layout(rect=[0, 0, 1, 0.96])  # Adjust layout to fit the title
    plt.show()

def plot_velocity(t, Vcggs):
    
    titles = ['East', 'North', 'Zenith', 'Total Speed']
    ylimits = [[0, 500], [-150, 0], [-200, 800],[0,1000]]
    yticks = [
        np.arange(0, 501, 50),
        np.arange(-150, 1, 15),
        np.arange(-200, 801, 100),
        np.arange(0, 1001, 100),
    ]
    
    # Debugging- Validate input dimensions
    #if Vcggs.shape[1] != 3:
    #    raise ValueError("Vcggs must have exactly 3 columns (x, y, z velocities).")

    # Create 2x2 subplots
    fig, ax = plt.subplots(2, 2, sharex=True)
    ax = ax.flatten()  # Flatten to access subplots as a 1D array
    
    # Plot each coordinate position over time
    for i in range(3):
        ax[i].plot(t, Vcggs[:, i], 'k', linestyle='None', markersize=8, marker='.')
        ax[i].set_title(titles[i], fontsize=20, fontweight='bold', fontname='Arial')
        ax[i].set_ylabel('Speed (m/s)',fontsize=16, fontweight='bold', fontname='Arial')
        ax[i].set_ylim(ylimits[i])
        ax[i].set_yticks(yticks[i])
        ax[i].grid()
        
    # Plot total displacement (norm)
    total_displacement = np.linalg.norm(Vcggs, axis=1)
    ax[3].plot(t, total_displacement, 'k', linestyle='None', markersize=8, marker='.')
    ax[3].set_title(titles[3], fontsize=20, fontweight='bold', fontname='Arial')
    ax[3].set_ylabel('Speed (m/s)',fontsize=16, fontweight='bold', fontname='Arial')
    ax[3].set_ylim(ylimits[3])
    ax[3].set_yticks(yticks[3])
    ax[3].grid()
    
    # Set x-axis label on the bottom plots
    ax[2].set_xlabel('Time (s)',fontsize=16, fontweight='bold', fontname='Arial')
    ax[3].set_xlabel('Time (s)',fontsize=16, fontweight='bold', fontname='Arial')
    ax[2].set_xlim([0, 70])
    ax[2].set_xticks(np.linspace(0, 70, 7))
    

    fig.tight_layout(rect=[0, 0, 1, 0.96])  # Adjust layout to fit the title
    plt.show()

def plot_altitude(t, Rcggs, C):
    n= len(t)
    x= np.zeros((3,n))

    #Altitude, linear range, and haversine range
    for k in range(n):
        Rcge = C['GS']['Rgse'] + Rcggs[k, :]  # Absolute vehicle position WRT Earth in ENZ coordinates
        rcge = np.linalg.norm(Rcge)           # Vehicle range WRT Earth
        rgse = np.linalg.norm(C['GS']['Rgse'])  # Ground station range WRT Earth

        x[0, k] = rcge - C['Re']  # Altitude above mean equator
        x[1, k] = np.linalg.norm(Rcggs[k, :])  # Linear range
        theta = np.arccos(np.dot(Rcge, C['GS']['Rgse']) / (rcge * rgse))  # Angular displacement
        x[2, k] = theta * C['Re']  # Haversine range
    
    plt.figure()
    plt.plot(t, x[0, :], 'k', linestyle='None', markersize=8, marker='.')
    plt.title('Altitude Above Mean Equator', fontsize=20, fontweight='bold', fontname='Arial')
    plt.xlabel('Time (s)', fontsize=16, fontweight='bold', fontname='Arial')
    plt.ylabel('Altitude (m)', fontsize=16, fontweight='bold', fontname='Arial')
    plt.xlim([0, 70])
    plt.ylim([0, 6000])
    plt.xticks(np.arange(0, 71, 7), fontsize=12, fontname='Arial')
    plt.yticks(np.arange(0, 6001, 600), fontsize=12, fontname='Arial')
    plt.grid(True)
    plt.show()

def plot_acceleration(t,S,Rcggs,Vcggs,C):
    n= len(t)
    Acge= np.zeros((3,n))
    # []Allocates memory for the vehicle's inertial accelerations.

    a= np.zeros(n)

    for k in range(n):
        dSdt = projectile_eom(t[k], S[:, k], C)
        # [m/s, m/s^2] Evaluates the equations of motion.

        Acggs = dSdt[3:6]
        # [m/s^2] Vehicle acceleration relative to the ground station in ENZ coordinates.

        We = np.asarray(C["GS"]["We"]).reshape(3)  # Earth's angular velocity vector
        Rcggs_k = np.asarray(Rcggs[k, :]).reshape(3)  # Position vector
        Vcggs_k = np.asarray(Vcggs[k, :]).reshape(3)  # Velocity vector

        Acge[:, k] = (
        C["GS"]["Agse"]
        + np.cross(We, np.cross(We, Rcggs_k))
        + 2 * np.cross(We, Vcggs_k)
        + Acggs
        ) / C["g"]
        # [m/s^2] Inertial vehicle acceleration WRT the Earth in ENZ coordinates.

        a[k] = np.linalg.norm(Acge[:, k])
        # [g's] Inertial vehicle acceleration WRT the Earth.

    plt.figure()
    plt.plot(t, a, 'k', linestyle='None', markersize=8, marker='.')
    plt.title('Inertial Acceleration', fontsize=20, fontweight='bold', fontname='Arial')
    plt.xlabel('Time (s)', fontsize=16, fontweight='bold', fontname='Arial')
    plt.ylabel('Acceleration (g\'s)', fontsize=16, fontweight='bold', fontname='Arial')
    plt.xlim([0, 70])
    plt.ylim([0, 20])
    plt.xticks(np.arange(0, 71, 7), fontsize=12, fontname='Arial')
    plt.yticks(np.arange(0, 21, 2), fontsize=12, fontname='Arial')
    plt.grid(True)
    plt.show()

def plot_angles(t,S):
    n = len(t)
    azimuth = np.zeros(n)  # Initialize 1D array to store azimuth
    elevation = np.zeros(n)  # Initialize 1D array to store elevation

    for k in range(n):
        Rcggs = S[0:3, k]  # Extract the first three elements of the k-th column
        e = Rcggs[0]  # East component
        n = Rcggs[1]  # North component
        z = Rcggs[2]  # Zenith component

        # Calculate azimuth angle in degrees
        azimuth[k] = np.arctan2(e, n) * 180 / np.pi

        # Calculate the magnitude of Rcggs (avoid division by zero)
        norm_Rcggs = np.linalg.norm(Rcggs)

        # Avoid invalid values for arcsin (clip the value between -1 and 1)
        if norm_Rcggs != 0:  # Prevent division by zero
            elevation[k] = np.arcsin(np.clip(z / norm_Rcggs, -1, 1)) * 180 / np.pi
        else:
            elevation[k] = 0
    
    # Create figure
    fig, axs = plt.subplots(2, 1, sharex=True, tight_layout=True)

    # Titles and limits
    Titles = ['Azimuth', 'Elevation']
    YLimits = [[-180, 180], [-90, 90]]
    YTicks = [np.arange(-180, 181, 36), np.arange(-90, 91, 18)]

    for k in range(2):
        axs[k].plot(t, azimuth if k == 0 else elevation, 'k.', markersize=8, linewidth=1)
        axs[k].set_title(Titles[k], fontsize=16, fontweight='bold', fontname='Arial')
        axs[k].set_xlabel('Time (s)', fontsize=12, fontweight='bold', fontname='Arial')
        axs[k].set_ylabel(Titles[k], fontsize=12, fontweight='bold', fontname='Arial')
        axs[k].set_ylim(YLimits[k])
        axs[k].set_yticks(YTicks[k])
        axs[k].set_xlim([0, 70])
        axs[k].grid(True)

    plt.show()
    

def plot_pressure(t, S, C):
    n = len(t)
    dynamic_pressure = np.zeros(n)  # Initialize array for dynamic pressure

    for k in range(n):
        Rcggs = S[0:3, k]  # Position vector at time step k
        Vcggs = S[3:6, k]  # Velocity vector at time step k

        # Recalculate Rcge, Vcge, and Vatm at each step to find Vinf
        Rcge = C['GS']['Rgse'] + Rcggs  # Absolute vehicle position WRT the Earth
        Vcge = C['GS']['Vgse'] + np.cross(C['GS']['We'], Rcggs) + Vcggs  # Absolute vehicle velocity
        Vatm = np.cross(C['GS']['We'], Rcge)  # Atmosphere velocity
        Vinf = Vcge - Vatm  # True air velocity

        # Get density using altitude (magnitude of Rcggs as altitude)
        _, _, D = standard_atmosphere(np.linalg.norm(Rcggs))

        # Calculate the magnitude of true air velocity
        Vinf_mag = np.linalg.norm(Vinf)

        # Calculate dynamic pressure [kPa]
        dynamic_pressure[k] = 0.5 * D * Vinf_mag**2 / 1000  # Convert to kPa

    # Plot dynamic pressure vs. time
    plt.figure(figsize=(10, 6), facecolor='w')
    plt.plot(t, dynamic_pressure, 'k.', markersize=8, label="Dynamic Pressure")
    plt.title('Dynamic Pressure', fontsize=20, fontweight='bold', fontname='Arial')
    plt.xlabel('Time (s)', fontsize=16, fontweight='bold', fontname='Arial')
    plt.ylabel('Pressure (kPa)', fontsize=16, fontweight='bold', fontname='Arial')
    plt.xlim([0, 70])  
    plt.ylim([0, 500]) 
    plt.xticks(np.arange(0, 71, 7), fontsize=12, fontweight='bold', fontname='Arial')
    plt.yticks(fontsize=12, fontweight='bold', fontname='Arial')
    plt.grid(True) 
    plt.show()

def run_simulation():

    start_time = time.time()
    # Start time of the simulation 

    # Define the constants and initial conditions
    C = constants()

    # Initial state vector: position and velocity (east, north, zenith)
    initial_position = C['CG']['Rcggs']
    initial_velocity = C['CG']['Vcggs']
    S0 = np.hstack((initial_position, initial_velocity))

    # Define the time span for the simulation (0 to 70 seconds)
    t_span = (0, 70)
    t_eval = np.linspace(t_span[0], t_span[1], 500)  # Evaluate at 500 points within the time span

    # Run the simulation using solve_ivp
    solution = solve_ivp(
        fun=lambda t, S: projectile_eom(t, S, C),
        t_span=t_span,
        y0=S0,
        t_eval=t_eval,
        )

    # Extract the time and position data from the solution
    t = solution.t
    S = solution.y  
    Rcggs = solution.y[:3, :].T  # Extract position (first 3 components of S)
    Vcggs = solution.y[3:, :].T  # Last 3 components are velocity
    
    # Plot the results
    plot_position(t, Rcggs)
    plot_velocity(t,Vcggs)
    plot_altitude(t,Rcggs,C)
    plot_acceleration(t,S,Rcggs,Vcggs,C)
    plot_angles(t,S)
    plot_pressure(t,S,C)

    end_time = time.time()
    # End time of the simulation
    
    print("Simulation Complete!", round(end_time-start_time,3), "seconds") 

# Run the simulation
run_simulation()