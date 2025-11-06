#! /bin/bash

source activate pdal

python ../run_model.py \
    -i ./mimics/mimics_inputs_2024-02-26_Angola_P09.csv \
    -c ../model/data/ \
    -s ../model/code/