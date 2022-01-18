# run_bedtools_coverage.sh
# run as ./run_bedtools_coverage.sh cell1 C CAST dir_bowtie_result dir_bed

#!/bin/sh

# # cell1, cell2, ...
# file_name=$1
# # C,S
# flag=$2
# # CAST,129S1
# bed=$3
# # dir_in="/home/ygli/gam_paper_data/gam_seq_mapped/bowtie_result/"${file_name}"/"
# # dir_bed="/home/ygli/gam_paper_data/gam_seq_mapped/"${bed}".bed"

# dir_in=$4${file_name}"/"
# dir_bed=$5${bed}".bed"

# cd ${dir_in}
# for each_np in *_${flag}.sam
# do
#     np_name=${each_np%%'.sam'}
#     echo $np_name
#     samtools view -f 0x2 -bS ${np_name}.sam > ${np_name}.bam
#     samtools sort ${np_name}.bam > ${np_name}.sort.bam
#     samtools index ${np_name}.sort.bam
#     java -jar /home/ygli/tools/picard.jar MarkDuplicates REMOVE_DUPLICATES=True I=${np_name}.sort.bam O=${np_name}.sort.remove_dump.bam M=${np_name}.sort.marked_dup_metrics.txt
#     bedtools coverage -a ${dir_bed} -b ${np_name}.sort.remove_dump.bam > ${np_name}.depth.txt
# done


# ./run_bedtools_coverage.sh Cell7 mm10 mm10 "/home/ygli/gam_paper_data/gam_seq_mapped2/bowtie_result/" "/home/ygli/gam_paper_data/gam_seq_mapped2/" 
# for re calculate in 2021/4/29
# cell1, cell2, ...
file_name=$1
# C,S
flag=$2 #mm10
# CAST,129S1
bed=$3
# dir_in="/home/ygli/gam_paper_data/gam_seq_mapped/bowtie_result/"${file_name}"/"
# dir_bed="/home/ygli/gam_paper_data/gam_seq_mapped/"${bed}".bed"

dir_in=$4${file_name}"/"
dir_bed=$5${bed}".bed"

cd ${dir_in}
# for each_np in *_${flag}.sort.remove_dump.bam
for each_np in *_${flag}.sam
do
    np_name=${each_np%%'.sam'}
    echo $np_name
    samtools view -f 0x2 -bS ${np_name}.sam > ${np_name}.bam
    samtools sort ${np_name}.bam > ${np_name}.sort.bam
    samtools index ${np_name}.sort.bam
    java -jar /home/ygli/tools/picard.jar MarkDuplicates REMOVE_DUPLICATES=True I=${np_name}.sort.bam O=${np_name}.sort.remove_dump.bam M=${np_name}.sort.marked_dup_metrics.txt
    bedtools coverage -a ${dir_bed} -b ${np_name}.sort.remove_dump.bam > ${np_name}.depth.txt
done
