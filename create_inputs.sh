for i in $(seq 1 2 10)
do 
	jd=`echo $i| awk '{printf("%g", $1*1e11)}'`
	mkdir j_dens$i
	cd j_dens$i
cat<<EOF >input
# Material parameters:
Ms=477e3
L=200e-9
Ly=50e-9
Lz=5e-9
Aex=1.05e-11
alpha=0.02

Nwires=1


# Field's parameters:
H=0 #A/m
f=500e6

# Simulation's parameters:
Dt=1e-13
totaltime=10e-9
out_time=0.1e-9
integrator=RK4

#Pinning field:
pin=0
Temperature=0
#spin torque paramters
P=0.7
beta=0.04
j=$jd #A/m^2

program=Benchmark2
#program=Bifurcation
EOF
cp ../EXEC . 
./EXEC
cd ..
done 
