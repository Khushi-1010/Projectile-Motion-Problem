# Projectile-Motion-Problem
MAE 3405 Flight Dynamics Honors Python Project

##Projectile-Motion Simulator (Python Edition)

You’re looking at a full-Python rewrite of the MAE 3405 “Mighty Mouse” projectile-motion homework. While the standard assignment prescribed MATLAB, this repo delivers the same deliverables—plus a few extras—in an open, reproducible Python package.

## ✅ Original Problem Statement (From Homework #6)

At a coordinated universal time of **UTC = 2024:02:28:21:00:00**, the students of MAE 3405 launch a prototype version of the **“Mighty Mouse” projectile** from the middle of the **50-yard line of AT&T Stadium**:

- Ground station latitude: **ϕgs = 32° 44’ 52’’ N**
- Ground station longitude: **λgs = 97° 5’ 34’’ W**
- Ground station altitude: **hgs = 185 m**

The projectile launch conditions (defined **with respect to the launch site**) are:

- Launch speed: **827 m/s**
- Azimuth: **107°**
- Elevation: **60°**

### What the assignment asks
Simulate the projectile’s trajectory in **ENZ coordinates** (East–North–Zenith) with respect to the ground station while considering:

- **Gravity**
- **Atmospheric drag**

The simulation runs for a time period equal to the time required for the projectile to **impact the initial altitude above mean equator** after launch (i.e., stop when the projectile returns to its initial altitude level).

The projectile is modeled as a sphere with:

- Mass: **mcg = 35 kg**
- Drag coefficient: **cd = 0.82**
- Profile area: **ar = (6.00225 × 10⁻³)π m²**

##Why This Exists:
- Honors credit ➜ deeper dive – I rebuilt the solver from scratch to explore numerical integrators, Earth-rotation corrections, and a modular plotting pipeline.
- MATLAB → Python – No license walls; everything runs with nothing more exotic than NumPy/SciPy/Matplotlib.
- Teaching aid – The code is heavily commented and split into bite-sized functions so classmates (and future students) can compare MATLAB vs Python workflows.

##Results You’ll Get:
- Altitude Above Mean Equator vs Time – obeys the spec’s fonts & tick spacing.
- 2×2 Displacement suite – East, North, Zenith, and ground-range.
- 2×2 Velocity suite – same breakdown for speeds.
- Inertial Acceleration (g’s) – reveals max g-loading ~12 g.
- Azimuth & Elevation tracks – antenna-pointing angles.
- Dynamic Pressure curve – bonus insight.

##License & Attribution
© 2024 Khushi Piparava. For academic or commercial use, please reach out first (krp3505@mavs.uta.edu).
