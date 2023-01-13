#!/bin/bash

python get_data.py
./clustering.sh
python blast_model.py