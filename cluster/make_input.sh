#!/bin/sh
# Generates a DWPC input file for one realisation of the 3-wire Ni laser study.
# usage: make_input.sh <base_input> <seed> <laser|ref> > input
# All physics parameters are set here; edit them in one place.
BASE=$1; SEED=$2; MODE=$3
if [ "$MODE" = laser ]; then USE_TTM=1; else USE_TTM=0; fi

TOTALTIME=1.2e-9      # s
T_PULSE=0.2e-9        # s, time of the pulse maximum
OUT_TIME=1e-12        # s
DT=5e-17              # s
P0=2.3279e22          # W/m^3, peak absorbed power (1.24e9 J/m^3 in a 50 fs FWHM pulse)
SIGMA=2.1233e-14      # s (50 fs FWHM)
TAU_M=200e-15         # s
TAU_COOL=0            # s (0 = no cooling)
RIJ=${RIJ:-50e-9}     # m, edge-to-edge gap between wires (env var); RIJ=inf -> no coupling
if [ "$RIJ" = inf ]; then COUPLING=0; RIJ=50e-9; else COUPLING=1; fi

sed -e 's/^Ms=.*/Ms=509407.45/' -e 's/^Aex=.*/Aex=1.19752e-11/' \
    -e 's/^Nwires=.*/Nwires=3/' -e "s/^use_DW_coupling=.*/use_DW_coupling=$COUPLING/" -e "s/^rij=.*/rij=$RIJ/" -e "s/^rng_seed=.*/rng_seed=$SEED/" \
    -e 's/^Temperature=.*/Temperature=300/' -e 's/^use_Ms_T=.*/use_Ms_T=1/' -e 's/^Tc=.*/Tc=635/' \
    -e 's/^Ms_exponent=.*/Ms_exponent=0.341/' -e 's/^Ms_rescaling_exponent=.*/Ms_rescaling_exponent=2.322/' \
    -e 's/^Aex_exponent=.*/Aex_exponent=2/' -e 's/^pinning_exponent=.*/pinning_exponent=2/' \
    -e "s/^use_TTM=.*/use_TTM=$USE_TTM/" -e 's/^TTM_gamma_e=.*/TTM_gamma_e=6e3/' -e 's/^TTM_C_l=.*/TTM_C_l=3.9e6/' \
    -e 's/^TTM_G=.*/TTM_G=8e17/' -e 's/^TTM_fluence=.*/TTM_fluence=0/' \
    -e "s/^TTM_t0=.*/TTM_t0=$T_PULSE\nTTM_P0=$P0/" -e "s/^TTM_sigma=.*/TTM_sigma=$SIGMA/" \
    -e "s/^TTM_tau_m=.*/TTM_tau_m=$TAU_M/" -e "s/^TTM_tau_cool=.*/TTM_tau_cool=$TAU_COOL/" \
    -e "s/^Dt=.*/Dt=$DT/" -e "s/^totaltime=.*/totaltime=$TOTALTIME/" -e "s/^out_time=.*/out_time=$OUT_TIME/" \
    -e 's/^integrator=.*/integrator=RK4/' -e 's/^program=.*/program=Benchmark3/' "$BASE"
