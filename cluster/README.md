# DWPC laser ensemble on the Slurm cluster

Paired runs: every seed is simulated with and without the laser pulse, using the
same thermal-noise sequence (rng_seed), so dx = x_laser - x_ref isolates the
effect of the pulse on each wire.

1. On node 104, in ~/DWPC (branch with patches 0001-0007 applied):
       bash cluster/deploy.sh          # copies + builds DWPC on all nodes
2. Submit (100 seeds; each task runs 2 x ~40 s single-core simulations):
       sbatch --array=1-100 cluster/dwpc_ensemble.sbatch
3. When squeue is empty, gather and analyse on 104:
       bash cluster/collect.sh
       python3 cluster/analyze.py ~/dwpc_ens_all
   -> dx_summary.csv, dx_final.csv, dwpc_ensemble.png

All physics parameters are in cluster/make_input.sh. Seeds are 1000 + array index,
so new batches (e.g. --array=101-300) add new seeds without overwriting old ones.
Node list for deploy/collect: NODES="ip1 ip2 ..." bash cluster/deploy.sh
