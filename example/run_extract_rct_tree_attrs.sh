#!/bin/bash

source activate pdal

python ../extract_rct_tree_attrs.py \
    -m /data/TLS2/angola/2024-02-26_P09.RiSCAN/matrix \
    -t /data/TLS2/angola/2024-02-26_P09.RiSCAN/rct_tile50m_overlap10m/tile_index.dat \
    -r /data/TLS2/angola/2024-02-26_P09.RiSCAN/rct_tile50m_overlap10m/tiled \
    -c Angola \
    -p P09 \
    -d 2024-02-26 \
    -o ./ \
    -f ./