# merge_overlap.sh

#!/bin/sh
# cell_dir="/home/ygli/gam_paper_data/gam_seq_mapped2/overlap_result_chr1_10kb/"
# merge_output_dir="/home/ygli/gam_paper_data/gam_seq_mapped2/overlap_result_chr1_10kb/"

cell_list=$1
code_dir=$2
cell_dir=$3
merge_output_dir=$4
first_cell=$5
final_cell=$6

echo "doing merge in cell..."
for cell in `cat $cell_list`
do
    echo $cell
    python ${code_dir}combine_cover.py $cell S ${merge_output_dir} >${cell_dir}${cell}"/cover_S.log"
    python ${code_dir}combine_cover.py $cell C ${merge_output_dir} >${cell_dir}${cell}"/cover_C.log"
done

cd ${merge_output_dir}
mkdir combine
merge_output_dir=${merge_output_dir}"combine/"

echo "doing merge between cell..."

# for the first 3 column
cut -f 1,2,3 ${cell_dir}${first_cell}/${first_cell}_C_total_cover.txt > ${merge_output_dir}first3.txt
# merge first cell
awk -F'\t' '{print $NF}' ${cell_dir}${first_cell}/${first_cell}_C_total_cover.txt | paste ${merge_output_dir}first3.txt - >${merge_output_dir}${first_cell}_C_combine.txt

for cell in `cat $cell_list`
do
    if [ "$cell" != "$final_cell" ];then
        awk -F'\t' '{print $NF}' ${cell_dir}${cell:0:4}`expr ${cell:4:5} + 1`/${cell:0:4}`expr ${cell:4:5} + 1`_C_total_cover.txt | paste ${merge_output_dir}${cell}_C_combine.txt - >${merge_output_dir}${cell:0:4}`expr ${cell:4:5} + 1`_C_combine.txt
    fi
    # # merge Cell8
    # awk -F'\t' '{print $NF}' ${cell_dir}Cell8/Cell8_C_total_cover.txt | paste ${merge_output_dir}Cell7_C_combine.txt - >${merge_output_dir}Cell8_C_combine.txt
    # # merge Cell9
    # awk -F'\t' '{print $NF}' ${cell_dir}Cell9/Cell9_C_total_cover.txt | paste ${merge_output_dir}Cell8_C_combine.txt - >${merge_output_dir}Cell9_C_combine.txt
    # # merge Cell10
    # awk -F'\t' '{print $NF}' ${cell_dir}Cell10/Cell10_C_total_cover.txt | paste ${merge_output_dir}Cell9_C_combine.txt - >${merge_output_dir}Cell10_C_combine.txt
    # # merge Cell11
    # awk -F'\t' '{print $NF}' ${cell_dir}Cell11/Cell11_C_total_cover.txt | paste ${merge_output_dir}Cell10_C_combine.txt - >${merge_output_dir}Cell11_C_combine.txt
done

# for the first 3 column
cut -f 1,2,3 ${cell_dir}${first_cell}/${first_cell}_S_total_cover.txt > ${merge_output_dir}first3.txt
# merge first cell
awk -F'\t' '{print $NF}' ${cell_dir}${first_cell}/${first_cell}_S_total_cover.txt | paste ${merge_output_dir}first3.txt - >${merge_output_dir}${first_cell}_S_combine.txt

for cell in `cat $cell_list`
do
    if [ "$cell" != "$final_cell" ];then
        awk -F'\t' '{print $NF}' ${cell_dir}${cell:0:4}`expr ${cell:4:5} + 1`/${cell:0:4}`expr ${cell:4:5} + 1`_S_total_cover.txt | paste ${merge_output_dir}$cell_S_combine.txt - >${merge_output_dir}${cell:0:4}`expr ${cell:4:5} + 1`_S_combine.txt
    fi
done


# # for the first 3 column
# cut -f 1,2,3 ${cell_dir}Cell7/Cell7_S_total_cover.txt > ${merge_output_dir}first3.txt
# # merge Cell7
# awk -F'\t' '{print $NF}' ${cell_dir}Cell7/Cell7_S_total_cover.txt | paste ${merge_output_dir}first3.txt - >${merge_output_dir}Cell7_S_combine.txt
# # merge Cell8
# awk -F'\t' '{print $NF}' ${cell_dir}Cell8/Cell8_S_total_cover.txt | paste ${merge_output_dir}Cell7_S_combine.txt - >${merge_output_dir}Cell8_S_combine.txt
# # merge Cell9
# awk -F'\t' '{print $NF}' ${cell_dir}Cell9/Cell9_S_total_cover.txt | paste ${merge_output_dir}Cell8_S_combine.txt - >${merge_output_dir}Cell9_S_combine.txt
# # merge Cell10
# awk -F'\t' '{print $NF}' ${cell_dir}Cell10/Cell10_S_total_cover.txt | paste ${merge_output_dir}Cell9_S_combine.txt - >${merge_output_dir}Cell10_S_combine.txt
# # merge Cell11
# awk -F'\t' '{print $NF}' ${cell_dir}Cell11/Cell11_S_total_cover.txt | paste ${merge_output_dir}Cell10_S_combine.txt - >${merge_output_dir}Cell11_S_combine.txt
