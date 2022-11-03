#!/bin/sh
# run_sam2seg.sh
# conda activate bwa_env
# bash run_sam2seg.sh Cell1 mm10 CAST dir_bowtie_result dir_bed dir_hickit dir_SNP
# bash run_sam2seg.sh Cell7 mm10 mm10 "/home/ygli/gam_paper_data/gam_seq_mapped2/bowtie_result/" "/home/ygli/gam_paper_data/gam_seq_mapped_control/" "/home/ygli/tools/hickit-0.1_x64-linux" "/data01/ygli/SNP_hui/phased_SNP.tsv"
# for calculate snp of each reads

file_name=$1
flag=$2 #mm10
bed=$3
# dir_in="/home/ygli/gam_paper_data/gam_seq_mapped/bowtie_result/"${file_name}"/"
# dir_bed="/home/ygli/gam_paper_data/gam_seq_mapped/"${bed}".bed"

dir_in=$4${file_name}"/"
dir_bed=$5${bed}".bed"

dir_hickit=$6 
dir_SNP=$7

cd ${dir_in}
for each_np in *_${flag}.sort.remove_dump.bam
do
    np_name=${each_np%%'.sort.remove_dump.bam'}
    echo $np_name
    samtools view ${each_np} | sort > ${each_np}.sam.sort
    ${dir_hickit}/k8 ${dir_hickit}/hickit.js sam2seg -q 0 -d 0 -v ${dir_SNP} ${each_np}.sam.sort | ${dir_hickit}/k8 ${dir_hickit}/hickit.js chronly -y - > ${each_np}.sam.sort.callsnp
done
