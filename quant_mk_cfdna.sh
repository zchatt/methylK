#!/bin/bash
# This script assigns paired and trimmed sequencing reads from the same samples of tissue type "identify" to tissue type "interest" and "contrast" using the k-mer index produced by $methylK_dir/quant_mk_cell.sh.
# The C and G content of each PE read uniquely assigned is then counted.

##################
## Housekeeping ##
##################
# Check if necessary software is installed
type parallel >/dev/null 2>&1 || { echo >&2 "parallel is required but not installed.  Aborting."; exit 1; }
type samtools >/dev/null 2>&1 || { echo >&2 "samtools is required but not installed.  Aborting."; exit 1; }
type kallisto >/dev/null 2>&1 || { echo >&2 "kallisto is required but not installed.  Aborting."; exit 1; }
type trimmomatic >/dev/null 2>&1 || { echo >&2 "trimmomatic is required but not installed.  Aborting."; exit 1; }


# Check that all values have been supplied
if [[ $# -eq 0 ]] ; then
    echo 'No arguments supplied. Aborting.'
    exit 1
fi
if [ "$#" -ne 5 ]; then
  echo "The 5 required arguments not supplied. Aborting." >&2
  exit 1
fi

# Assign variables
if [ -z "$PS1" ]; then
    methylK_dir=$1 # methylK directory
    targets=$2 # phenotype dataframe
    sdir=$3 # location of .fastq files
    odir=$4 # output directory
    read_length=$5

fi

export njobs=4 #set this to a number so that minerva is happy

echo ""
echo "[ INPUTS ]"
echo ""
echo "methylK directory == $1"
echo "targets file == $2"
echo "sample directory (.fastq.gz) == $3"
echo "output directory == $4"
echo "read length == $5"

# ensure targets file properly formatted
perl -pi -e 's/\r\n/\n/g' $targets

#####################
### Read Trimming ###
#####################
cd $odir
# Ensure .fastq files are present and make lists for sample types "identify"
echo ""
echo "[ Checking all .fastq.gz files are present ]"
echo ""

rm read1 read2 2> /dev/null
for i in $(awk '$3 ~ /identify/ {print $1}' $targets); do
    # make read list
    ls $sdir/${i}.R1.fastq.gz >> read1
    ls $sdir/${i}.R2.fastq.gz >> read2
if [ ! -f $(ls -d $sdir/${i}.R1.fastq.gz) ]
then
    echo "$0: File '${i}' not found." >&2
  exit 1
fi
done

# create softlink to Nextera.fa
ln -sf $methylK_dir/NexteraPE-PE.fa NexteraPE-PE.fa

# run trimmomatic in parallel
parallel --xapply -j $njobs --eta trimmomatic PE {1} {2} {1/.}_truncated {1/.}_untruncated {2/.}_truncated {2/.}_untruncated ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25 CROP:$read_length :::: read1 :::: read2

# rename truncated reads
rename R1.fastq_paired.fq_truncated paired_truncated_R1.fastq *R1.fastq_paired.fq_truncated
rename R2.fastq_paired.fq_truncated paired_truncated_R2.fastq *R2.fastq_paired.fq_truncated
# gzip files
gzip *truncated_R1.fastq
gzip *truncated_R2.fastq

#######################################################
### Truncate paired and trimmed reads - deprecated  ###
#######################################################
# cd $odir
# # Ensure .fastq files are present and make lists for sample types "interest" & "contrast"
# echo ""
# echo "[ Checking all paired/trimmed .fastq.gz files are present ]"
# echo ""

# rm p_read1 p_read2 2> /dev/null
# for i in $(awk '$3 ~ /identify/ {print $1}' $targets); do
#     # make read list
#     ls $odir/${i}.R1.fastq_paired.fq.gz >> p_read1
#     ls $odir/${i}.R2.fastq_paired.fq.gz >> p_read2
# if [ ! -f $(ls -d $odir/${i}.R1.fastq_paired.fq.gz) ]
# then
#     echo "$0: File '${i}' not found." >&2
#   exit 1
# fi
# done

# # Truncate reads
# echo ""
# echo "[ Truncating .fastq files ]"
# echo ""

# # read length to truncate to
# R1_length=$read_length
# R2_length=$read_length

# # truncate reads
# #parallel --xapply -j+0 --eta fastqutils truncate {1} $R1_length '>' {1.}_truncated :::: p_read1
# #parallel --xapply -j+0 --eta fastqutils truncate {1} $R2_length '>' {1.}_truncated :::: p_read2
# # run trimmomatic in parallel
# parallel --xapply -j $njobs --eta trimmomatic CROP $read_length PE {1} {2} {1.}_truncated {1.}_untruncated {2.}_truncated {2.}_untruncated :::: p_read1 :::: p_read2

# # rename truncated reads
# rename R1.fastq_paired.fq_truncated paired_truncated_R1.fastq *R1.fastq_paired.fq_truncated
# rename R2.fastq_paired.fq_truncated paired_truncated_R2.fastq *R2.fastq_paired.fq_truncated
# # gzip files
# gzip *truncated_R1.fastq
# gzip *truncated_R2.fastq


############################################
### assign (pseudoalign) truncated reads ###
############################################
echo ""
echo "[ Assigning (pseudoalign) truncated reads  ]"
echo ""

# make directories and sample names to deposit results and make list of truncated paired-read R1/R2 fastq files
rm samp mdir pt_read1 pt_read2 2> /dev/null
cp /dev/null mdir
cp /dev/null samp
for i in $(awk '$3 ~ /identify/ {print $1}' $targets); do
    echo $odir/${i}"_out" >> mdir
    echo ${i} >> samp
    ls -d ${i}.paired_truncated_R1.fastq.gz >> pt_read1
    ls -d ${i}.paired_truncated_R2.fastq.gz >> pt_read2

if [ ! -f $(ls -d $odir/${i}.paired_truncated_R1.fastq.gz) ]
then
    echo "$0: File '${i}' not found." >&2
  exit 1
fi
done

# run kallisto quant in parallel
parallel --xapply -j $njobs --eta kallisto quant -i $odir/master_methylotype.kidx -o {3} --pseudobam {1} {2} :::: pt_read1 :::: pt_read2 :::: mdir

# convert pseudo.bam to pseudo.sam
for file in $(cat mdir); do
samtools view $file/pseudoalignments.bam >> $(basename ${file%_out}).pseudoalignments.sam
done

# count unique pseduoaligned reads and total pseudoaligned reads
echo ""
echo "[ Summarising assigned reads ]"
echo ""

for SAM in *pseudoalignments.sam; do
OUTNAME1=${SAM%%.pseudoalignments.sam}.psdccnt
OUTNAME2=${SAM%%.pseudoalignments.sam}_total.psdccnt
# select pseudoaligned reads from .sam file output
awk -F":f:|\t" '{if($3!="*" && $1 !~ /^@HD|^@PG|^@SQ/)print;}' $SAM > $odir/temp.tmp3
# remove duplicated reads
awk '!seen[$1]++' $odir/temp.tmp3 > $odir/tmp.uniq
# count uniquely pseudoaligned reads to one cells methylotpye from .sam file output
rm $odir/$OUTNAME1 2> /dev/null
rm $odir/$OUTNAME2 2> /dev/null
awk '$12 == "NH:i:1" {print $3}' $odir/tmp.uniq | sort | uniq -c >> $odir/$OUTNAME1
awk '{print $3}' $odir/tmp.uniq | awk '{split($1,a,":methylotype"); print a[1]}' | sort | uniq -c >> $odir/$OUTNAME2
done

##################################################
### count cg content in pseudoaligned PE reads ###
##################################################
# count methylated loci in unique pseudoaligned reads
cd $odir
awk '$3 ~ /identify/ {print $1}' $targets > sample_name_list
awk -F'[/]' '{print $(NF-0) "_out"}' sample_name_list > outdir_list
awk '$3 ~ /identify/ {print $1".pseudoalignments.sam"}' $targets  > psam_list
parallel --xapply -j 4 --eta $methylK_dir/upsd_cgcount.sh {1} {2} {3} :::: sample_name_list :::: outdir_list :::: psam_list 2> /dev/null

# move files to output directory
mv *out/*cgcount $odir

# cleanup tmp files
rm sample_name_list outdir_list psam_list




