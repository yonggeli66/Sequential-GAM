# run_all_chromatin_calculate_mapped_chr_process.sh
#!/bin/sh

dir_file="/home/ygli/gam_paper_data/gam_seq_mapped2/"

echo "doing calculate_mapped_chr_process.sh for all chromatin ... "
for chr in {chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chrX}
do
    echo $chr
    # bash /home/ygli/gam_paper_data/gam_seq_mapped2/calculate_mapped_chr_process.sh chr.list ${chr} 10
    bash calculate_mapped_chr_process.sh ${dir_file}/cell.list ${chr} 10
done

for chr in {chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chrX}
do
    echo $chr
    echo "Cell11_C_overlap"
    cat "/home/ygli/gam_paper_data/gam_seq_mapped2/overlap_result_"$chr"_10kb/combine/Cell11_C_overlap.table"
    echo "Cell11_S_overlap"
    cat "/home/ygli/gam_paper_data/gam_seq_mapped2/overlap_result_"$chr"_10kb/combine/Cell11_S_overlap.table"
done