#!/bin/bash

set -euo pipefail

input_path=""
output_path="output"
file_type=".pdf"

while getopts :i:o:t:h opt
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
    t)
      file_type=$OPTARG
      ;;
    h)
      cat <<EOU
-i path to input file
-o path to output file
-t type output file will be
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

if [[  $input_path == "" ]]; then
   echo "Need input file! (see -h)"
   false
fi

image=/lustre/scratch117/cellgen/cellgeni/TIC-misc/tic-1129/actions/images/cluster_investigation_0_1.sif
mount_options="/lustre,/nfs"

/software/singularity-v3.5.3/bin/singularity exec -B $mount_options $image \
img2pdf $input_path -o $output_path$file_type
