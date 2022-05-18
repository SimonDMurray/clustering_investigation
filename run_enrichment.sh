#!/bin/bash

set -euo pipefail

im=/lustre/scratch117/cellgen/cellgeni/TIC-misc/tic-1129/actions/images/cluster_investigation_0_1.sif

script=$1
seurat=$2
rcl=$3
hallmark=$4
cells=$5
path=$6

/software/singularity-v3.5.3/bin/singularity exec -B /lustre -B /nfs $im Rscript $script --seurat=$seurat --rcl=$rcl --hallmark=$hallmark --cellnumber=$cells --output=$path
