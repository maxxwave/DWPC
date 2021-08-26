CC=g++
CFLAGS= -std=c++11 -O3 #-Wall

SRCFILES= src/main.cpp \
	  src/initialize.cpp \
	  src/calculate.cpp \
	  src/storage_variables.cpp \
	  src/runge_kutta4th.cpp \
	  src/euler_integrator.cpp \
	  src/bifurcation.cpp \
	  src/benchmark.cpp\
	  src/rc.cpp
EXE=DWPC

${EXE} : ${SRCFILES}
	${CC} $^ ${CFLAGS}  -o ${EXE} -llapack -g -pg
