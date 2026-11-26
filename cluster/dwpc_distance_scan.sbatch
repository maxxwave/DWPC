#!/bin/bash
# Distance scan: array index -> (gap, seed). With NGAP gaps and NSEED seeds:
#   sbatch --array=0-$((NGAP*NSEED-1)) cluster/dwpc_distance_scan.sbatch
# e.g. 8 gaps x 100 seeds:  sbatch --array=0-799 cluster/dwpc_distance_scan.sbatch
# Results: ~/dwpc_scan/gap_<rij>/run_<seed>/{laser,ref}/  (gather with collect.sh SRC=dwpc_scan)
#SBATCH --job-name=dwpc_scan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=dwpc_scan_%A_%a.log
set -e
GAPS=(25e-9 50e-9 100e-9 150e-9 250e-9 400e-9 800e-9 inf)   # edge-to-edge gap; inf = coupling off
NSEED=${NSEED:-100}
SEED_BASE=1000
IG=$((SLURM_ARRAY_TASK_ID / NSEED)); IS=$((SLURM_ARRAY_TASK_ID % NSEED))
RIJ=${GAPS[$IG]}; SEED=$((SEED_BASE + IS))
DWPC_DIR=$HOME/DWPC
RUN_DIR=$HOME/dwpc_scan/gap_$RIJ/run_$SEED
echo "task $SLURM_ARRAY_TASK_ID gap $RIJ seed $SEED on $(hostname)"
for MODE in laser ref; do
    mkdir -p "$RUN_DIR/$MODE"; cd "$RUN_DIR/$MODE"
    RIJ=$RIJ sh "$DWPC_DIR/cluster/make_input.sh" "$DWPC_DIR/input" "$SEED" "$MODE" > input
    "$DWPC_DIR/DWPC" > log 2> err || true
    if [ ! -s Benchmark3.data ] || grep -q "^Error" err; then echo "run failed"; tail -5 err; exit 1; fi
done
hostname > "$RUN_DIR/host"
