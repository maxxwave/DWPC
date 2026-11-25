#!/bin/bash
# Run on node 104: gathers ~/$SRC (default dwpc_ens) from all nodes into ~/${SRC}_all
set -e
NODES="${NODES:-172.29.250.102 172.29.250.105 172.29.250.106 172.29.250.107 172.29.250.108}"
SRC=${SRC:-dwpc_ens}; DST=~/${SRC}_all
mkdir -p $DST
rsync -a ~/$SRC/ $DST/ 2>/dev/null || true
for N in $NODES; do
    rsync -a -e "ssh -p 2222 -o BatchMode=yes" uaic@$N:$SRC/ $DST/ || echo "no results on $N"
done
find $DST -name Benchmark3.data | wc -l
