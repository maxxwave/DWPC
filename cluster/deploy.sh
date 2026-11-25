#!/bin/bash
# Run on node 104 as uaic, from ~/DWPC (patches applied, cluster/ folder present).
# Copies the DWPC sources to every node over ssh -p 2222 and builds them there.
set -e
NODES="${NODES:-172.29.250.102 172.29.250.105 172.29.250.106 172.29.250.107 172.29.250.108}"
SSH="ssh -p 2222 -o BatchMode=yes"

# build locally first
sudo apt-get install -y -qq g++ liblapack-dev > /dev/null
make -C ~/DWPC clean >/dev/null 2>&1 || true
make -C ~/DWPC

for N in $NODES; do
    echo "=== $N"
    rsync -az --delete -e "$SSH" --exclude .git --exclude DWPC --exclude '*.o' ~/DWPC/ uaic@$N:DWPC/
    $SSH uaic@$N 'sudo apt-get install -y -qq g++ liblapack-dev >/dev/null && make -C ~/DWPC >/dev/null && ls -la ~/DWPC/DWPC'
done
