#!/bin/bash 
# The following scripts were used within "Chatterton et al, Methods for detecting Brain-Cell derived Cell-Free DNA, 2020" for the preparation of bisulfite genome (genome_prepare.sh) alignment and methylation calling using bismark (meth_trim_align_call.sh).

################## 
## Housekeeping ##
##################
# Check if necessary software is installed
type trimmomatic >/dev/null 2>&1 || { echo >&2 "trimmomatic is required but not installed.  Aborting."; exit 1; }
type parallel >/dev/null 2>&1 || { echo >&2 "parallel is required but not installed.  Aborting."; exit 1; }
type bismark >/dev/null 2>&1 || { echo >&2 "bismark is required but not installed.  Aborting."; exit 1; }
type samtools >/dev/null 2>&1 || { echo >&2 "samtools is required but not installed.  Aborting."; exit 1; }
hash $PICARD_JAR >/dev/null 2>&1 || { echo >&2 "picard is required but not installed.  Aborting."; exit 1; }
type bowtie2 >/dev/null 2>&1 || { echo >&2 "bowtie2 is required but not installed.  Aborting."; exit 1; }

# Check that all values have been supplied
if [[ $# -eq 0 ]] ; then
    echo 'No arguments supplied. Aborting.'
    exit 1
fi
if [ "$#" -lt 5 ]; then
  echo "The 5 required arguments not supplied. Aborting." >&2
  exit 1
fi

# Assign variables
if [ -z "$PS1" ]; then
	methylK_dir=$1 # location of methylK github repository
	genome=$2 # genome .fa file
	sdir=$3 # sample directory containing all PE .fastq files 
	odir=$4 # ouput directory
	targets=$5 # phenotype dataframe
fi

# path for bowtie2, found automatically 
bowtie2=$(dirname $(which bowtie2))

#set this to a number so that the server is happy
export njobs=4

if [ $6 == "identify" ]
  then
  	echo 'Running on sample types "identify"'
		# ensure R1 and R2 .fastq files are present and make read lists for sample types "identify"
		rm read1
		for FILES in $(awk '$3 ~ /interest|contrast/ {print $1}' $targets)
		do
		i=$sdir/${FILES}.R1.fastq.gz
		ls -d ${i} >> read1
		if [ ! -f $(ls -d $sdir/${FILES}.R1.fastq.gz) ]
		then
		echo "$0: File '${i}' not found." >&2
		exit 1
		fi
		done

		rm read2
		for FILES in $(awk '$3 ~ /interest|contrast/ {print $1}' $targets)
		do
		i=$sdir/${FILES}.R2.fastq.gz
		ls -d ${i} >> read2
		if [ ! -f $(ls -d $sdir/${FILES}.R2.fastq.gz) ]
		then
		echo "$0: File '${i}' not found." >&2
		exit 1
		fi
		done

		###########
		### RUN ###
		###########
		# running from the output folder 
		# Prepare bisulfite genome from .fasta, index and create a dictionary necessary for .vcf to .fasta conversion
		mkdir ./bismark_genome
		ln -sf $genome ./bismark_genome/
		$methylK_dir/genome_prepare.sh ./bismark_genome/

		# Perform bisulfite sequencing alignment and DNA methylation calling
		$methylK_dir/meth_trim_align_call.sh ./bismark_genome/ $bowtie2 $sdir $methylK_dir

		# cleanup
		rm bam_list out_dir output_log

	else
		 echo 'Running on sample types "interest" & "contrast"'
		# ensure R1 and R2 .fastq files are present and make read lists for sample types "interest" & "contrast"
		rm read1
		for FILES in $(awk '$3 ~ /interest|contrast/ {print $1}' $targets)
		do
		i=$sdir/${FILES}.R1.fastq.gz
		ls -d ${i} >> read1
		if [ ! -f $(ls -d $sdir/${FILES}.R1.fastq.gz) ]
		then
		echo "$0: File '${i}' not found." >&2
		exit 1
		fi
		done

		rm read2
		for FILES in $(awk '$3 ~ /interest|contrast/ {print $1}' $targets)
		do
		i=$sdir/${FILES}.R2.fastq.gz
		ls -d ${i} >> read2
		if [ ! -f $(ls -d $sdir/${FILES}.R2.fastq.gz) ]
		then
		echo "$0: File '${i}' not found." >&2
		exit 1
		fi
		done

		###########
		### RUN ###
		###########
		# running from the output folder 
		# Prepare bisulfite genome from .fasta, index and create a dictionary necessary for .vcf to .fasta conversion
		mkdir ./bismark_genome
		ln -sf $genome ./bismark_genome/
		$methylK_dir/genome_prepare.sh ./bismark_genome/

		# Perform bisulfite sequencing alignment and DNA methylation calling
		$methylK_dir/meth_trim_align_call.sh ./bismark_genome/ $bowtie2 $sdir $methylK_dir

		# cleanup
		rm bam_list out_dir output_log

fi
