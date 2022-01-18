# merge_9_cell.sh

#!/bin/sh
# eg. ./merge_9_cell.sh chr1 10 Cell11 Cell20


# cell1_dir="/home/ygli/gam_paper_data/gam_seq_mapped2/overlap_result_chr1_10kb/"
# cell2_dir="/home/ygli/gam_paper_data/gam_seq_mapped_219/overlap_result_chr1_10kb/"

chr=$1
resolution=$2
cell1_name=$3
cell2_name=$4

cell1_dir="/home/ygli/gam_paper_data/gam_seq_mapped2/overlap_result_"$chr"_"$resolution"kb/combine/"
cell2_dir="/home/ygli/gam_paper_data/gam_seq_mapped_219/overlap_result_"$chr"_"$resolution"kb/combine/"
merge_output_dir="/home/ygli/gam_paper_data/gam_seq_mapped_219/merge_result_2_219/"

# cd merge_output_dir
# mkdir $output_file_name

flag="C"
file1_name=${cell1_name}"_"${flag}"_combine.txt"
file2_name=${cell2_name}"_"${flag}"_combine.txt"
output_file_name=$chr"_"$resolution"kb_"${flag}"_combine.txt"
awk -F'\t' '{for(i=1;i<=3;i++){$i=""}; print $0,'\t'}' ${cell2_dir}${file2_name} | paste ${cell1_dir}${file1_name} - >${merge_output_dir}${output_file_name}

flag="S"
file1_name=${cell1_name}"_"${flag}"_combine.txt"
file2_name=${cell2_name}"_"${flag}"_combine.txt"
output_file_name=$chr"_"$resolution"kb_"${flag}"_combine.txt"
awk -F'\t' '{for(i=1;i<=3;i++){$i=""}; print $0,'\t'}' ${cell2_dir}${file2_name} | paste ${cell1_dir}${file1_name} - >${merge_output_dir}${output_file_name}

