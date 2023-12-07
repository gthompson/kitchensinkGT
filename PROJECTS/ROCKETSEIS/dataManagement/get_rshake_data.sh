#!/bin/bash
DATADIR="$HOME/DATA/KSC_RSHAKE"
echo $DATADIR
#rsync -avzr -e "ssh -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null" --progress myshake@raspberryshake.local:/opt/data ./
rsync -av myshake@128.217.252.92:/opt/data/archive/ $DATADIR
rsync -av myshake@128.217.250.104:/opt/data/archive/ $DATADIR
