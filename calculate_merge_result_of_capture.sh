# calculate_merge_result_of_capture.sh

code_dir="/home/ygli/gam_paper_data/code/"

cd $code_dir
bash run_merge_9_cell_mm10.sh

dir_file="/home/ygli/gam_paper_data/gam_seq_mapped_219/merge_result_2_219/"
bash run_calculate_calculate_sparse.sh ${code_dir} ${dir_file}

bash run_convert_combine_to_data_direct_for_combine_mm10.sh

cd dir_file
makir set_1cell_1 set_1cell_2 set_1cell_3 set_1cell_4 set_1cell_5 set_1cell_6 set_1cell_7 set_1cell_8 set_1cell_9
for cell_list in {set_1cell_1,set_1cell_2,set_1cell_3,set_1cell_4,set_1cell_5,set_1cell_6,set_1cell_7,set_1cell_8,set_1cell_9}
do
    bash run_convert_combine_to_data_direct_for_set_mm10.sh ${cell_list}".list"
done
