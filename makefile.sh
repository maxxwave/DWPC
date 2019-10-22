#!/bin/bash

g++ -O3   src/main.cpp src/initialize.cpp src/calculate.cpp src/storage_variables.cpp src/euler_integrator.cpp src/runge_kutta4th.cpp src/bifurcation.cpp -o EXEC
