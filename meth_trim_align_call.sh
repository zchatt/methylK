#!/bin/bash
# The following scripts perform trimming, alignment and methylation calling from .fastq files

# reference genome directory
BISGENOME=$1

# path to bowtie 2
bowtie2=$2

#path to the folder that contains the fq files
sdir=$3

#path to methylK
methylK_dir=$4

# create softlink to Nextera.fa
ln -sf $methylK_dir/NexteraPE-PE.fa NexteraPE-PE.fa

#####################
### Read Trimming ###
#####################

# run trimmomatic in parallel
parallel --xapply -j $njobs --eta trimmomatic PE {1} {2} {1/.}_paired.fq.gz {1/.}_unpaired.fq.gz {2/.}_paired.fq.gz {2/.}_unpaired.fq.gz ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25 :::: read1 :::: read2

#################
### Alignment ###
#################
# make read 1 list
ls -d *R1*_paired* > p_read1

# make read 2 list
ls -d *R2*_paired* > p_read2

# make output directories for each sample
awk -F'[/]' '{print $NF "_out"}' p_read1 > out_dir

# run bismark in parallel
parallel --xapply -j $njobs --eta bismark --bowtie2 --non_directional --bam --path_to_bowtie $bowtie2 -o {3} --temp_dir {3} $BISGENOME -1 {1} -2 {2} :::: p_read1 :::: p_read2 :::: out_dir > output_log

######################################
### bismark methylation exctractor ###
######################################
# make .bam file list
ls -d */*bam > bam_list

# run bismark methylation extraction in parallel
mkdir -p results
parallel --xapply -j $njobs --eta bismark_methylation_extractor -p --comprehensive --merge_non_CpG --cytosine_report --CX -o results --genome_folder $BISGENOME {1} :::: bam_list

# cleanup files
rm read1 read2 p_read1 p_read2

# cleanup directory
mkdir -p bismark_results
mv *out/* bismark_results
mv results/* bismark_results
rm -r results *out
