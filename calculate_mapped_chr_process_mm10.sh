# calculate_mapped_chr_process_mm10.sh

#!/bin/sh
# ./calculate_mapped_chr_process.sh cell_list chr resolution(bp)
# eg1
# bash calculate_mapped_chr_process_mm10.sh /home/ygli/gam_paper_data/gam_seq_mapped_219/cell.list chr1 10
# eg2
# for resolution in {1000,10000,100000,1000000}
# do
#     for chr_each in {1..19}
#     do
#         echo ${resolution}_chr${chr_each}
#         bash calculate_mapped_chr_process_mm10.sh /home/ygli/gam_paper_data/gam_seq_mapped_219/cell.list chr${chr_each} ${resolution}
#     done
# done

cell_list=$1
chromatin=$2
resolution=$3

# code_dir="/home/ygli/gam_paper_data/code_upload/"
# dir_gam_seq_raw="/home/ygli/gam_paper_data/gam_seq_raw_219/"
# dir_data_after_adaptor="/home/ygli/gam_paper_data/gam_seq_mapped_219/data_after_adaptor/"
# dir_bowtie_result="/home/ygli/gam_paper_data/gam_seq_mapped_219/bowtie_result/"
# dir_data="/home/ygli/gam_paper_data/gam_seq_mapped_219/"
# dir_bed="/home/ygli/gam_paper_data/gam_seq_mapped_219/"
# first_cell="Cell17"
# final_cell="Cell20"

code_dir="/home/ygli/gam_paper_data/code_upload/"
dir_gam_seq_raw="/home/ygli/gam_paper_data/gam_seq_raw_414/"
dir_data_after_adaptor="/home/ygli/gam_paper_data/gam_seq_mapped_414/data_after_adaptor/"
dir_bowtie_result="/home/ygli/gam_paper_data/gam_seq_mapped_414/bowtie_result/"
dir_data="/home/ygli/gam_paper_data/gam_seq_mapped_414/"
dir_bed="/home/ygli/gam_paper_data/gam_seq_mapped_414/"
first_cell="Cell25"
final_cell="Cell39"

# code_dir="/home/ygli/gam_paper_data/code/"
# dir_gam_seq_raw="/home/ygli/gam_paper_data/gam_seq_raw2/"
# dir_data_after_adaptor="/home/ygli/gam_paper_data/gam_seq_mapped2/data_after_adaptor/"
# dir_bowtie_result="/home/ygli/gam_paper_data/gam_seq_mapped2/bowtie_result/"
# dir_data="/home/ygli/gam_paper_data/gam_seq_mapped2/"
# dir_bed="/home/ygli/gam_paper_data/gam_seq_mapped2/"
# first_cell="Cell7"
# final_cell="Cell11"

echo "calculate the map of each cell in "$cell_list" over "$chromatin" in "$resolution" bp"

echo "overlap each cell ..."
cd $code_dir
mkdir -p ${dir_data}overlap_result_"$chromatin"_"$resolution"
for each_cell in `cat $cell_list`
do
    echo $each_cell
    bash overlap_np.sh $dir_bowtie_result $each_cell "mm10" $dir_bowtie_result"chrom_windows/mm10_"$chromatin"_windows_"$resolution".bed" $dir_data"overlap_result_"$chromatin"_"$resolution"/" 
done

# merge the result of overlap
echo "merge the result of each cell ..."
cd $code_dir
bash merge_overlap_mm10.sh $cell_list $code_dir $dir_data"overlap_result_"$chromatin"_"$resolution"/" $dir_data"overlap_result_"$chromatin"_"$resolution"/" $first_cell $final_cell

# calculate the sparsity
echo "calculate the sparsity ... "
cd $code_dir
python calculate_sparse.py ${final_cell}"_mm10_combine.txt" $dir_data"overlap_result_"$chromatin"_"$resolution"/combine/"

cd $dir_data"overlap_result_"$chromatin"_"$resolution"/combine/"
echo "calculate the map of each cell in "$cell_list" over "$chromatin" in "$resolution
echo ${final_cell}"_mm10_overlap.table"
cat ${final_cell}_mm10_combine_count.txt | awk -F'\t' '{print $NF}'| sort | uniq -c |tee ${final_cell}_mm10_overlap.table


