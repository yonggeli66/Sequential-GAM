#!/bin/sh
# run_calculate_calculate_sparse.sh

code_dir=$1
dir_file=$2

cd $dir_file
for file in *_combine.txt
do
    echo ${file}
    python ${code_dir}calculate_sparse.py ${file} ${dir_file}
done