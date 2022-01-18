# calculate_mapped_chr_process.sh

#!/bin/sh
# ./calculate_mapped_chr_process.sh cell_list chr resolution(kb)
# eg ./calculate_mapped_chr_process.sh /home/ygli/gam_paper_data/gam_seq_mapped_219/cell.list chr1 10

cell_list=$1
chromatin=$2
resolution=$3

code_dir="/home/ygli/gam_paper_data/gam_seq_mapped_219/"
dir_gam_seq_raw="/home/ygli/gam_paper_data/gam_seq_raw_219/"
dir_data_after_adaptor="/home/ygli/gam_paper_data/gam_seq_mapped_219/data_after_adaptor/" 
dir_bowtie_result="/home/ygli/gam_paper_data/gam_seq_mapped_219/bowtie_result/" 
dir_bed="/home/ygli/gam_paper_data/gam_seq_mapped_219/"
first_cell="Cell17"
final_cell="Cell20"

echo "calculate the map of each cell in "$cell_list" over "$chromatin" in "$resolution" kb"

echo "overlap each cell ..."
cd $code_dir
mkdir overlap_result_"$chromatin"_"$resolution"kb
for each_cell in `cat $cell_list`
do
    echo $each_cell
    ./overlap_np.sh $dir_bowtie_result $each_cell S $dir_bowtie_result"chrom_windows/129S1_"$chromatin"_windows_"$resolution"kb.bed" $code_dir"overlap_result_"$chromatin"_"$resolution"kb/" 
    ./overlap_np.sh $dir_bowtie_result $each_cell C $dir_bowtie_result"chrom_windows/CAST_"$chromatin"_windows_"$resolution"kb.bed" $code_dir"overlap_result_"$chromatin"_"$resolution"kb/"
done

# merge the result of overlap
echo "merge the result of each cell ..."
cd $code_dir
./merge_overlap.sh $cell_list $code_dir $code_dir"overlap_result_"$chromatin"_"$resolution"kb/" $code_dir"overlap_result_"$chromatin"_"$resolution"kb/" $first_cell $final_cell

# calculate the sparsity
echo "calculate the sparsity ... "
cd $code_dir

python calculate_sparse.py ${final_cell}"_C_combine.txt" $code_dir"overlap_result_"$chromatin"_"$resolution"kb/combine/"
python calculate_sparse.py ${final_cell}"_S_combine.txt" $code_dir"overlap_result_"$chromatin"_"$resolution"kb/combine/"

cd $code_dir"overlap_result_"$chromatin"_"$resolution"kb/combine/"

echo "calculate the map of each cell in "$cell_list" over "$chromatin" in "$resolution" kb"
echo ${final_cell}"_C_overlap.table"
cat ${final_cell}_C_combine_count.txt | awk -F'\t' '{print $NF}'| sort | uniq -c |tee ${final_cell}_C_overlap.table
echo ${final_cell}"_S_overlap.table"
cat ${final_cell}_S_combine_count.txt | awk -F'\t' '{print $NF}'| sort | uniq -c |tee ${final_cell}_S_overlap.table
