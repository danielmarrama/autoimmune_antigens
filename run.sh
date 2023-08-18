#!/bin/bash

python pull_data.py
./clustering.sh
python blast_model.py