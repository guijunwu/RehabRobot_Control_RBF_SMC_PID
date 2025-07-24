# RBF, SMC, and PID Control for 2D Upper-Limb Rehabilitation Robot

This repository contains the MATLAB simulation code and demonstration video for the paper titled:

**"Robust Tracking Control in Rehabilitation Robotics: ST-SMC, RBF, and PID"**

## 📁 Folder Structure

- `code/`: MATLAB source code for the three controllers (PID, SMC, RBF).
- `video/`: Demonstration videos showing simulation or hardware experiments.

## 📌 Description

This study compares three control methods for a 2D upper-limb rehabilitation robot under external disturbances:

- **PID Control** (with adaptive gain)
- **Super-Twisting Sliding Mode Control (ST-SMC)**
- **RBF Neural Network Adaptive Control**

The simulation includes friction, sinusoidal disturbances, and admittance target models. All methods are evaluated using MATLAB simulation.

## ▶️ How to Run

1.Open the desired simulation script in the code/ folder using MATLAB.
2.Set the controller selection flags (e.g., PID, SMC, RBF) as needed.
3.Run the script to execute the simulation and generate results.
4.Use the provided plotting functions or scripts to visualize velocity tracking error, motor force, and other performance indicators.

## 📂 File Overview

| File Name                                                       | Description                                                                                                                                     |
| ------------------------------------------------ | ------------------------------------------------------------------------------------------------------ |
| `controller_pid.m`                                        | Adaptive PID controller implementation for velocity tracking                                             |
| `controller_st_smc.m`                                  | Super-Twisting Sliding Mode Controller (ST-SMC) implementation                                  |
| `controller_rbf.m`                                         | Radial Basis Function Neural Network (RBF NN) adaptive controller implementation    |
| `evaluation_disturbance10.m`                   | Simulation and evaluation script under mild disturbance (amplitude = 10)                      |
| `evaluation_disturbance20.m`                   | Simulation and evaluation script under stronger disturbance (amplitude = 20)              |
| `combined_methods_disturbance10.m` | Comparative simulation of all three controllers under disturbance amplitude = 10       |
| `combined_methods_disturbance20.m` | Comparative simulation of all three controllers under disturbance amplitude = 20       |
| `combined_methods_with_tanh.m`         | Modified simulation using `tanh()` instead of `sign()` to suppress motor chattering        |
| `demo_hardware_test.mp4`                      | Experimental validation video demonstrating real system performance                           |
| `README.md`                                               | Project documentation and usage instructions (this file)                                                      |


## 📽️ Demonstration Videos

Experimental validation results on the physical system.

## 📎 License

For academic use only. If you use the code or video in your research, please cite our work.
