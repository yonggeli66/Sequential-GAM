# overlap_np.sh

#!/bin/sh

file_dir=$1
file_name=$2
flag=$3
chr_bed_dir=$4
dir_out=$5

cd ${dir_out}
mkdir ${file_name}

dir_out=${dir_out}${file_name}"/"

cd ${file_dir}${file_name}
for each_np in *_${flag}.sort.remove_dump.bam
do
    echo ${each_np}
    # trans .remove_dump.bam into bed
    bamToBed -i ${each_np} > ${each_np%%'.bam'}.bed
    # do overlap
    bedtools intersect -a ${chr_bed_dir} -b ${each_np%%'.bam'}.bed -c >${dir_out}${each_np%%'.sort.remove_dump.bam'}.cover.bed
done
