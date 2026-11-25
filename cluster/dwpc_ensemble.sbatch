#!/bin/bash
# Slurm job array: one array task = one seed = one paired run (with / without laser).
#   sbatch --array=1-100 dwpc_ensemble.sbatch
# Each node needs ~/DWPC (patched + built, see deploy.sh). Results stay on the
# node that ran the task in ~/dwpc_ens/run_<seed>/ and are gathered by collect.sh.
#SBATCH --job-name=dwpc_ens
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=dwpc_ens_%A_%a.log

set -e
SEED_BASE=1000
SEED=$((SEED_BASE + SLURM_ARRAY_TASK_ID))
DWPC_DIR=$HOME/DWPC
RUN_DIR=$HOME/dwpc_ens/run_$SEED

echo "task $SLURM_ARRAY_TASK_ID seed $SEED on $(hostname)"
for MODE in laser ref; do
    mkdir -p "$RUN_DIR/$MODE"
    cd "$RUN_DIR/$MODE"
    sh "$DWPC_DIR/cluster/make_input.sh" "$DWPC_DIR/input" "$SEED" "$MODE" > input
    "$DWPC_DIR/DWPC" > log 2> err || true      # DWPC returns 1 on success for Benchmark3
    if ! grep -q "^Setting up" err || [ ! -s Benchmark3.data ]; then
        echo "run $SEED $MODE failed"; tail -5 err; exit 1
    fi
    grep -q "^Error" err && { echo "run $SEED $MODE stopped with an error"; grep "^Error" err; exit 1; }
done
echo "$(hostname)" > "$RUN_DIR/host"
echo "done seed $SEED"
