# run_convert_combine_to_data_direct_for_combine.sh
#!/bin/sh

code_dir="/home/ygli/gam_paper_data/gam_seq_mapped_219/"
dir_input_root="/home/ygli/gam_paper_data/gam_seq_mapped_219/merge_result_2_219/"
frequency_output=9
csv_dir="/home/ygli/gam_paper_data/gam_seq_mapped_219/merge_result_2_219/combine_B8_B16_calculate_with_index.csv"
dir_output_root="/home/ygli/gam_paper_data/gam_seq_mapped_219/merge_result_2_219/"

cd $dir_output_root
S_dir="all_chr_10kb_S"
C_dir="all_chr_10kb_C"
mkdir $S_dir $C_dir

echo "convert combine to data that can be directly input to calculate... "
for chr in {chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chrX}
do
    echo $chr
    chr_name=${chr}"_10kb"
    python ${code_dir}convert_combine_to_data_direct_for_combine.py ${frequency_output} "S" ${dir_input_root} $csv_dir $chr_name ${dir_output_root}${S_dir}"/"
    python ${code_dir}convert_combine_to_data_direct_for_combine.py ${frequency_output} "C" ${dir_input_root} $csv_dir $chr_name ${dir_output_root}${C_dir}"/"
done