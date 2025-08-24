# MPC and MOAS Implementation for Aircraft Model

## Overview
This repository contains MATLAB scripts implementing:
1. Model Predictive Control (MPC) for a linearized aircraft model.
2. Approximation of Maximal Output Admissible Set (MOAS) for the system outputs.

The MPC is implemented manually without using toolboxes, including constraints on:
- Input (elevator angle)
- Rate of change of input
- Output (pitch angle and altitude)

MOAS is approximated by checking output constraint satisfaction over a finite prediction horizon.

## Files
- `MPC_Matrices.m` : Builds MPC cost function and constraint matrices.
- `moasApprox.m`   : Approximates the MOAS for given output constraints.
- `main.m`         : Example script demonstrating MPC simulation and MOAS visualization.

## Requirements
- MATLAB R2023a or newer.

## Purpose
Educational project demonstrating hands-on implementation of MPC and MOAS for control systems.
