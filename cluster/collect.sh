#!/bin/bash
# Run on node 104: gathers ~/dwpc_ens/run_* from all nodes into ~/dwpc_ens_all/
set -e
NODES="${NODES:-172.29.250.102 172.29.250.105 172.29.250.106 172.29.250.107 172.29.250.108}"
mkdir -p ~/dwpc_ens_all
rsync -a ~/dwpc_ens/ ~/dwpc_ens_all/ 2>/dev/null || true
for N in $NODES; do
    rsync -a -e "ssh -p 2222 -o BatchMode=yes" uaic@$N:dwpc_ens/ ~/dwpc_ens_all/ || echo "no results on $N"
done
ls -d ~/dwpc_ens_all/run_* | wc -l
