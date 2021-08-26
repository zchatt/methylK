#!/bin/bash
cd /project/DementiaCFDNA/Epigenetics/methylK_debug

git clone https://github.com/zchatt/methylK.git
chmod +x methylK/*

module load trimmomatic/0.36
module load parallel/20160222
module load bismark/0.18.2
module load samtools/1.9
module load picard/2.7.1
module load gatk/3.6.0 # required by vcf2fasta.sh
module load kallisto
module load python/2.7.9
module load R/3.6.3
module load bowtie2/2.2.5
module load bedtools
module load bwa
GATK_JAR=gatk
PICARD_JAR=picard
njobs=4

# Inputs for testing
# note - the methylK/test directory contains .fastq files and targets.txt that can be used to test scripts
methylK_dir=$(readlink -f methylK)
genome=$methylK_dir/tNGBS_n33_lambda1.3.fa # reference genome
sdir=$methylK_dir/test # directory containing all PE .fastq files 
targets=$methylK_dir/test/targets.txt # targets file
odir=$methylK_dir/test/output # ouput directory is where all results will be written

# Inputs for testing relative to above paths
mkdir $odir
snr_thresh_path=$odir # location of SNR threshold files. For convenience we also supply the pre-computed Glia and Neuron SNR thresholds within $methylK_dir/GN_thresh_snr
genome_bismark=$odir/bismark_genome/Bisulfite_Genome
bed_location=$odir/bismark_results

# hard-coded thresholds & values
cov_threshold=5 # coverage threshold in order to include cytosine in analysis
bin_threshold=50 # methylation % threshold for binarizing DNA methylation i.e. < 50% = T, > 50% = C
read_length=25 # length that reads will be truncated. This should be the shortest read length

# Run #
cd $odir

$methylK_dir/fastq_to_bed.sh $methylK_dir $genome $sdir $odir $targets

$methylK_dir/bismark_to_fasta.sh $odir $targets $cov_threshold $genome $genome_bismark $bed_location $bin_threshold $GATK_JAR $PICARD_JAR

$methylK_dir/quant_mk_cell.sh $methylK_dir $targets $sdir $odir $read_length

Rscript --vanilla $methylK_dir/snr_calc.R $odir $methylK_dir $targets

$methylK_dir/quant_mk_cfdna.sh $methylK_dir $targets $sdir $odir $read_length

Rscript --vanilla $methylK_dir/snr_quant.R $odir $methylK_dir $snr_thresh_path $targets
