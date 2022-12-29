# isi_proj_AGV

Realization and simulation of an AGV dynamic model for estimating the state space with Extended Kalman Filter and Unscented Kalman Filter.
The project has been realised for the course of Uncertain Systems Identification at University of Pisa (Course Of Studies in Robotics).

# Requirements
MATLAB 2022b or later versions.

# Instructions

Execute the script 'init.m' first to inizialize parameters.
Run the simulink model.
Run 'EKF.m' and 'UKF.m' to use the two filters for extimation.
Run 'plot_filters.m' to see the estimating result plots.
To see the regularized trajectory run 'EKF_for_smoother.m' and then 'RTS.m' to see plots.
