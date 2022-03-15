#!/bin/bash

iter=50
threads=8
tilesize=2048
PSF=PSF_dapi.tif

for file in dapi*tif; do
    [ -e "$file" ] || continue
    cmd="dw --iter $iter --tilesize $tilesize --threads $threads "$file" $PSF"
    echo "Will run: $cmd"
    eval $cmd
done
