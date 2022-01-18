# remove_adaptor.sh
#!/bin/sh
# remove_adaptor.sh file_name "/home/ygli/gam_paper_data/gam_seq_raw/" "/home/ygli/gam_paper_data/gam_seq_mapped/data_after_adaptor/"

file_name=$1
# dir_in="/home/ygli/gam_paper_data/gam_seq_raw/"${file_name}"/"
# dir_out="/home/ygli/gam_paper_data/gam_seq_mapped/data_after_adaptor/"${file_name}"/"

dir_in=$2${file_name}"/Cleandata/"
dir_out=$3${file_name}"/"

cd ${dir_in}"/"
for each_np in *
do
    fa1=${each_np}"/*R1.fq.gz"
    fast1=$(ls ${fa1})
    fa2=${each_np}"/*R2.fq.gz"
    fast2=$(ls ${fa2})
    trim_galore -q 20 --stringency 3 --length 20 --nextera -o ${dir_out} -e 0.1 --gzip --paired ${dir_in}"/"${fast1} ${dir_in}"/"${fast2}
done


