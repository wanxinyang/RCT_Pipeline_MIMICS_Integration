#! /bin/bash

source activate pdal

python ../create_parameters.py \
    -t 2024-02-26_Angola_P09_treeLevel_attr_rct_tile50m_overlap10m.csv \
    -r /data/TLS2/angola/2024-02-26_P09.RiSCAN/rct_tile50m_overlap10m/tiled \
    -c angola_P09_matched_stems.csv \
    -w wood_density.csv \
    -odir ./mimics/ \
    -a 1.05 \
    --min_branch_order 3 \
    --max_branch_order 5 \
    -d 10 \
    -H 3 \
    -f 0.43 1.2 5.4 \
    -s 0.25 0.5 0.75
