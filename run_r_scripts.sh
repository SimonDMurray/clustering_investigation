#!/bin/bash

set -euo pipefail

script=""
input_path=""
rds_type=false
output_path="./"
want_neighbours=false

while getopts :s:i:o:rnh opt
do
    case "$opt" in
    s)
      script=$OPTARG
      ;;
    i)
      input_path=$OPTARG
      ;;
    o)
      output_path=$OPTARG
      ;;
    r)
      rds_type=true
      ;;
    n)
      want_neighbours=true
      ;;
    h)
      cat <<EOU
-s script to run
-i path to input file
-r rds input type (tenx input type assumed without this flag)
-o path to output location (default CWD)
-n want neighbour outputs from FindNeighbours function
-h displays this message!
EOU
      exit
            ;;
    :) echo "Flag $OPTARG needs argument"
        exit 1;;
    ?) echo "Flag $OPTARG unknown"              
        exit 1;;
   esac
done

if [[  $script == "" ]]; then
   echo "Need script! (see -h)"
   false
fi

if [[  $input_path == "" ]]; then
   echo "Need input file! (see -h)"
   false
fi

if $rds_type; then
   echo "rds input mode selected"
   rds_var="-rds"
else
   rds_var=""
fi

if $want_neighbours; then
   echo "neighbour files requested"
   neighbours_var="-neighbours"
else
   neighbours_var="" 
fi

image=/lustre/scratch117/cellgen/cellgeni/TIC-misc/tic-1129/actions/images/cluster_investigation_v0.1.sif
mount_options="/lustre,/nfs"

/software/singularity-v3.5.3/bin/singularity exec -B $mount_options $image \
Rscript $script --input=$input_path --output=$output_path $rds_var $neighbours_var
