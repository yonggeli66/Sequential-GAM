# Sequential-GAM
Pipeline code for Sequential-GAM(Genome Architecture Mapping).

## mapping

`whole_preprocess.sh` include the whole processing of mapping.

 usage: ` bash whole_preprocess.sh Cell_name`

1.removing the adaptor by trim_galore

parameters are:

`  trim_galore -q 20 --stringency 3 --length 20 --nextera -o ${dir_out} -e 0.1 --gzip --paired ${dir_in}"Cleandata/"${fast1} ${dir_in}"Cleandata/"${fast2} `



2.mapping reads by bowtie2

parameters are:

`bowtie2 -p 10 -5 15 -x /home/ygli/index/CAST -1 B2F_NP1_R1.fq.gz -2 B2F_NP1_R2.fq.gz -S B2F_NP1_C.sam`



3.removing duplicates with samtools and Picard Toolkit

parameters are:

`samtools view -f 0x2 -bSB10F-NP252_S.sam > B10F-NP252_S.bam `

`samtools sort B10F-NP252_S.bam > B10F-NP252_S.sort.bam `

`samtools index B10F-NP252_S.sort.bam `

`java -jar /home/ygli/tools/picard.jar MarkDuplicates REMOVE_DUPLICATES=True I=B10F-NP252_C.sort.bam O=B10F-NP252_C.sort.remove_dump.bam M=B10F-NP252_C.sort.marked_dup_metrics.txt`



## capture table

1.calculating the capture or not for each window with given resolution (10kb in codes)

`bash run_all_chromatin_calculate_mapped_chr_process.sh`



2.merge capture situation in each cell into one table

`bash calculate_merge_result_of_capture.sh`



## model

The parameters of build_3D_structure_snpstrain_with_parameter.py are:

--chr_name: str,  default = "chr1", the chromosome to build 3D structure for.
--cell_name: str, default= "6", the cell name to build.
--method: str, default="overlap", the method to integrate the loci.
--my_step: int, default=180000, the step for optimizing to iteration.
--my_lr: float, default=1e-4, the learning rate for optimizing to search.
--strain_plot: str, default="C", the strain of genome.
--limit_LAD: float, default=1e-7.
--limit_res: float, default=1e+9.
--constant_LAD: float, default=1e-8.
--constant_res: float, default=1e-3.
--limit_power_LAD: int, default=3.
--limit_power_res: int, default=3.
--save_file_prefix: str,default="CS_depart_parameter_all_chr_10kb_", the file_prefix for save.

