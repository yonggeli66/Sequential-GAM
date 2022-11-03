# Sequential-GAM
Pipeline for Sequential-GAM(Genome Architecture Mapping).



## Contents

- [Overview](#overview)
- [Software Requirements](#Software Requirements)
- [mapping](#mapping)
- [merged table](#merged table)
- [run model](#run model)





## Overview



3D genome conformation plays an essential role in cell identity and gene regulation. However, how genome-wide hierarchical geometric structure is organized in single cells remains elusive. Here we developed Sequential-GAM (Sequential genome architecture mapping) to capture contiguous thin sections of the nucleus. Using Sequential-GAM, we constructed the hierarchical geometric structure and estimated the radial position of chromosomes, compartments, subcompartments, and genes in single cells. In this package, we provide the processing flow of the sequencing data and the construction of the 3D structure of individual chromosomes with the processed data.



## Software Requirements

Python 2.7

[trim_galore 0.6.4](https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md)

[Bowtie 2.3.5](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

[Picard tool](https://broadinstitute.github.io/picard/)

[dip-C](https://github.com/lh3/hickit)

[tensorflow 1.13.1](https://github.com/tensorflow/docs/tree/r1.13/site/en/api_docs)

[pandas](https://pandas.pydata.org/)

[scipy](https://scipy.org/)

[seaborn](https://seaborn.pydata.org/)



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



## merged table

1.calculating the capture or not for each window with given resolution (10kb in codes)

`bash run_all_chromatin_calculate_mapped_chr_process.sh`



2.merge capture situation in each cell into one table

`bash calculate_merge_result_of_capture.sh`



## run model

The parameters of build_3D_structure_snpstrain_with_parameter.py are:

--chr_name: str,  default = "chr1", the chromosome to build 3D structure for.
--cell_name: str, default= "6", the cell name to build.
--my_step: int, default=180000, the step for optimizing to iteration.
--my_lr: float, default=1e-4, the learning rate for optimizing to search.
--strain_plot: str, default="C", the strain of genome.
--limit_LAD: float, default=1e-3.
--limit_res: float, default=1e-3.
--constant_LAD: float, default=1e-8.
--constant_res: float, default=1e-3.
--limit_power_LAD: int, default=3.
--limit_power_res: int, default=1.
--save_file_prefix: str,default="newmm10_depart_parameter_all_chr_10kb_", the file_prefix for save.



Usage:

`python build_3D_structure_snpstrain_with_parameter_withchangedir.py --chr_name chr1 \
    --save_file_prefix newmm10_depart_parameter_all_chr_10kb_ \
    --cell_name 17 --strain_plot C \
    --cLAD_dropout 0.1 --limit_LAD 1e-7 --limit_res 1e+9 \
    --limit_power_LAD 3 --limit_power_res 3 \
    --power_law_p sampler \
    --input_z_file_prefix ./seqGAM_data/gam_seq_mapped/gam_seq_mapped_414_ \
--input_z_file_prefix ./seqGAM_data/gam_seq_mapped/gam_seq_mapped_414_ \
--input_LAD_file_prefix ./seqGAM_data/GSE17051_cLAD_regions_mm10_10000/ \
--powerlaw_paramter_file ./fit_parameter_each_chr_sampler.csv
`

Output:

The estimated 3D structure of each loci:

`/processed_data/set_1cell_28/newmm10_depart_parameter_all_chr_10kb_C/chr1_10kb_location_data_overlap.txt`


<!-- Notice:

CAST_snps_indels.tsv and 129S1_snps_indels.tsv in this github is incomplete. The complete file need to be download from ftp://ftp-mouse.sanger.ac.uk/REL-1807-SNPs_Indels/ -->



Expected run time:

The runtime in a computer with 264 GB RAM, 56 cores@2.60GHz for the demo is less than 2h.

