
for i in $(seq 1 1 5)
	do 
		#mkdir f$i
		cd f$i
		cp ../EXEC .
		c3=`echo $i|awk '{printf("%g", 0.638/$1)}'`
		c0=`echo $c3|awk '{printf("%e", 2.21117e-21/$1)}'`
		c1=`echo $c3|awk '{printf("%e", -3.8271e-15/$1)}'`
		c2=`echo $c3|awk '{printf("%e", -1.2866e-6/$1)}'`
		c4=`echo $c3|awk '{printf("%e", 1.632e8/$1)}'`

cat<<EOF >input

# Material parameters:
Ms=477e3
L=200e-9
Ly=50e-9
Lz=5e-9
Aex=1.05e-11
alpha=0.02

Nwires=1

a0=$c0
a1=$c1
a2=$c2
a3=$c3
a4=$c4
a5=0
a6=0
a7=0
a8=0
# Field's parameters:
H=800 #A/m
f=500e6

# Simulation's parameters:
Dt=1e-13
totaltime=20e-9
out_time=0.1e-9
integrator=RK4

#Pinning field:
pin=0
program=Bifurcation
EOF

./EXEC

cd ..
done
