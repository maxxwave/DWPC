CC=g++
CFLAGS=-O3 

SRCFILES= src/main.cpp \
	  src/initialize.cpp \
	  src/calculate.cpp \
	  src/storage_variables.cpp \
	  src/euler_integrator.cpp \
	  src/runge_kutta4th.cpp \
	  src/bifurcation.cpp \
	  src/rc.cpp
EXE=EXEC

${EXE} : ${SRCFILES}
	${CC} $^ ${CFLAGS}  -o ${EXE}
