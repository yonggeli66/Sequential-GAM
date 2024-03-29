# whole_preprocess.sh
# include the all preprocess, for the detail see the map_preprocess.md
# bash whole_preprocess.sh file_name
# e.g. bash whole_preprocess.sh Cell1
# e.g. bash whole_preprocess.sh Cell_control
#!/bin/sh

file_name=$1

# edit dir here
# code_dir="/home/ygli/gam_paper_data/gam_seq_mapped4/"
# fastqScreen_dir="/home/ygli/tools/fastq_screen_v0.14.0/"
# dir_gam_seq_raw="/home/ygli/gam_paper_data/gam_seq_raw4/"
# dir_bowtie_result="/home/ygli/gam_paper_data/gam_seq_mapped4/bowtie_result/" 
# dir_bed="/home/ygli/gam_paper_data/gam_seq_mapped4/"

# code_dir="/home/ygli/gam_paper_data/code/"
# fastqScreen_dir="/home/ygli/tools/fastq_screen_v0.14.0/"
# dir_gam_seq_raw="/home/ygli/gam_paper_data/gam_seq_raw2/"
# dir_data_after_adaptor="/home/ygli/gam_paper_data/gam_seq_mapped2/data_after_adaptor/" 
# dir_bowtie_result="/home/ygli/gam_paper_data/gam_seq_mapped2/bowtie_result/" 
# dir_bed="/home/ygli/gam_paper_data/gam_seq_mapped2/"

# code_dir="/home/ygli/gam_paper_data/code/"
# fastqScreen_dir="/home/ygli/tools/fastq_screen_v0.14.0/"
# dir_gam_seq_raw="/home/ygli/gam_paper_data/gam_seq_raw_219/"
# dir_data_after_adaptor="/home/ygli/gam_paper_data/gam_seq_mapped_219/data_after_adaptor/" 
# dir_bowtie_result="/home/ygli/gam_paper_data/gam_seq_mapped_219/bowtie_result/" 
# dir_bed="/home/ygli/gam_paper_data/gam_seq_mapped_219/"

# code_dir="/home/ygli/gam_paper_data/code/"
# fastqScreen_dir="/home/ygli/tools/fastq_screen_v0.14.0/"
# dir_gam_seq_raw="/home/ygli/gam_paper_data/gam_seq_raw_control/"
# dir_data_after_adaptor="/home/ygli/gam_paper_data/gam_seq_mapped_control/data_after_adaptor/" 
# dir_bowtie_result="/home/ygli/gam_paper_data/gam_seq_mapped_control/bowtie_result/" 
# dir_bed="/home/ygli/gam_paper_data/gam_seq_mapped_control/"

# code_dir="/home/ygli/gam_paper_data/code_upload/"
# fastqScreen_dir="/home/ygli/tools/fastq_screen_v0.14.0/"
# dir_gam_seq_raw="/home/ygli/gam_paper_data/gam_seq_raw_222/"
# dir_data_after_adaptor="/home/ygli/gam_paper_data/gam_seq_mapped_222/data_after_adaptor/" 
# dir_bowtie_result="/home/ygli/gam_paper_data/gam_seq_mapped_222/bowtie_result/" 
# dir_bed="/home/ygli/gam_paper_data/gam_seq_mapped_222/"

code_dir="/home/ygli/gam_paper_data/code_upload/"
fastqScreen_dir="/home/ygli/tools/fastq_screen_v0.14.0/"
dir_gam_seq_raw="/home/ygli/gam_paper_data/gam_seq_raw_414/"
dir_data_after_adaptor="/home/ygli/gam_paper_data/gam_seq_mapped_414/data_after_adaptor/" 
dir_bowtie_result="/home/ygli/gam_paper_data/gam_seq_mapped_414/bowtie_result/" 
dir_bed="/home/ygli/gam_paper_data/gam_seq_mapped_414/"

# code_dir="/home/ygli/gam_paper_data/code_upload/"
# fastqScreen_dir="/home/ygli/tools/fastq_screen_v0.14.0/"
# dir_gam_seq_raw="/home/ygli/gam_paper_data/gam_seq_raw_420/"
# dir_data_after_adaptor="/home/ygli/gam_paper_data/gam_seq_mapped_420/data_after_adaptor/" 
# dir_bowtie_result="/home/ygli/gam_paper_data/gam_seq_mapped_420/bowtie_result/" 
# dir_bed="/home/ygli/gam_paper_data/gam_seq_mapped_420/"


# code_dir="../"
# fastqScreen_dir="dir_to_fastq_screen/fastq_screen_v0.14.0/"
# dir_gam_seq_raw="dir_to_fastq_data/gam_seq_raw_control/"
# dir_data_after_adaptor="dir_to_save_data_removed_adaptor/data_after_adaptor/" 
# dir_bowtie_result="dir_to_save_mapped_data/bowtie_result/" 
# dir_bed="dir_to_reference_chr_size_with_bed_format/gam_seq_mapped_control/"

# remove adaptor
echo "=========remove adaptor========="
mkdir -p ${dir_data_after_adaptor}/${file_name}
cd ${code_dir}
# conda activate hic_env
bash remove_adaptor.sh ${file_name} ${dir_gam_seq_raw} ${dir_data_after_adaptor}

# bowtie mapping with script
echo "=========bowtie mapping with script========="
cd ${dir_bowtie_result}
mkdir ${file_name}
cd ${code_dir}
# ./run_bowtie.sh ${file_name} ${dir_data_after_adaptor} ${dir_bowtie_result} 2> ${dir_bowtie_result}"/"${file_name}"/mapping_result.txt"
bash run_bowtie_mm10.sh ${file_name} ${dir_data_after_adaptor} ${dir_bowtie_result} 2> ${dir_bowtie_result}"/"${file_name}"/"${file_name}"_mapping_result.txt"

# statistic
echo "=========statistic for coverage========="
# ./run_bedtools_coverage.sh ${file_name} C CAST ${dir_bowtie_result} ${dir_bed} 
# ./run_bedtools_coverage.sh ${file_name} S 129S1 ${dir_bowtie_result} ${dir_bed}
# python combine_depth.py ${file_name} C ${dir_bowtie_result}
# python combine_depth.py ${file_name} S ${dir_bowtie_result}
bash run_bedtools_coverage.sh ${file_name} mm10 mm10 ${dir_bowtie_result} ${dir_bed}
python combine_depth.py ${file_name} mm10 ${dir_bowtie_result}
# bash run_bedtools_coverage.sh Cell_control mm10 mm10 "/home/ygli/gam_paper_data/gam_seq_mapped_control/bowtie_result/" "/home/ygli/gam_paper_data/gam_seq_mapped_control/"

# calculate the fastqScreen
echo "=========fastqScreen========="
cd ${fastqScreen_dir}
mkdir ${file_name} 
cd ${file_name}
../fastq_screen --conf ../fastq_screen_my.conf ${dir_data_after_adaptor}${file_name}"/*.fq.gz"


