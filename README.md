# Projectile-Motion-Problem
MAE 3405 Flight Dynamics Honors Python Project

Projectile-Motion Simulator (Python Edition)

You’re looking at a full-Python rewrite of the MAE 3405 “Mighty Mouse” projectile-motion homework. While the standard assignment prescribed MATLAB, this repo delivers the same deliverables—plus a few extras—in an open, reproducible Python package.

Why This Exists:
- Honors credit ➜ deeper dive – I rebuilt the solver from scratch to explore numerical integrators, Earth-rotation corrections, and a modular plotting pipeline.
- MATLAB → Python – No license walls; everything runs with nothing more exotic than NumPy/SciPy/Matplotlib.
- Teaching aid – The code is heavily commented and split into bite-sized functions so classmates (and future students) can compare MATLAB vs Python workflows.

Results You’ll Get:
- Altitude Above Mean Equator vs Time – obeys the spec’s fonts & tick spacing.
- 2×2 Displacement suite – East, North, Zenith, and ground-range.
- 2×2 Velocity suite – same breakdown for speeds.
- Inertial Acceleration (g’s) – reveals max g-loading ~12 g.
- Azimuth & Elevation tracks – antenna-pointing angles.
- Dynamic Pressure curve – bonus insight.

License & Attribution
© 2024 Khushi Piparava. For academic or commercial use, please reach out first (krp3505 @ mavs.uta.edu).
