#!/bin/bash

for i in 1 5 10 15 20 25 50 
do
mkdir d$i
cd d$i

l=`echo $i | awk '{ printf("%g", $1*1e-9)}'`

cat<<EOF>input
# Material parameters:
Ms=477e3
L=1200e-9
Ly=50e-9
Lz=5e-9
Aex=1.05e-11
alpha=0.02

Nwires=1
#Energy profile
a0=0 #2.21117e-21
a1=0 #-3.8271e-15
a2=-3.2866e-6
a3=0 #0.61164
a4=2.632e8
a5=0
a6=0
a7=0
a8=0

# Field's parameters:
H=1500 #A/m
f=500e6

use_DW_coupling=0

#roughnss parameters
use_edge_roughness=1
H_edge_max=4000
edge_scale=$l
edge_seed=1234

# Simulation's parameters:
Dt=1e-13
totaltime=10e-9
out_time=0.1e-9
integrator=RK4

#Pinning field:
pin=1
Temperature=0
#spin torque paramters
P=0.7
beta=0.04
j=0 #1e12 #A/m^2
#program=multi_dw
program=Benchmark2
#program=RC
#program=spin_current
#program=Bifurcation
EOF

cp ../DWPC . 
./DWPC 
cd .. 
done
