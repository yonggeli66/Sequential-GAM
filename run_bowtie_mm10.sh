# run_bowtie_mm10.sh
#!/bin/sh
# run_bowtie.sh file_name "/home/ygli/gam_paper_data/gam_seq_mapped/data_after_adaptor/" "/home/ygli/gam_paper_data/gam_seq_mapped/bowtie_result/" 

file_name=$1
# dir_in="/home/ygli/gam_paper_data/gam_seq_mapped/data_after_adaptor/"${file_name}"/"
# dir_out="/home/ygli/gam_paper_data/gam_seq_mapped/bowtie_result/"${file_name}"/"

dir_in=$2${file_name}"/"
dir_out=$3${file_name}"/"

cd $dir_in

for each_np in *_R1_val_1.fq.gz
do
    echo
    np_name=${each_np%%'_R1_val_1.fq.gz'}
    echo $np_name" mapping mm10"
    bowtie2 -p 10 -5 15 -x /home/ygli/index/mm10 -1 ${np_name}_R1_val_1.fq.gz -2 ${np_name}_R2_val_2.fq.gz -S ${dir_out}${np_name}_mm10.sam
done

