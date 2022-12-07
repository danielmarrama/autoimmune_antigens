#!/bin/bash

cd cd-hit
./cd-hit -i ../../human_proteome.fasta -o ../human_db40 -c 0.4 -n 2 -d 0 -T 4
