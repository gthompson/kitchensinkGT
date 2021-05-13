#!/bin/bash
rsync -avzr -e "ssh -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null" --progress myshake@raspberryshake.local:/opt/data ./

