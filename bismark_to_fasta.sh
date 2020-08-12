#!/bin/bash
# This script converts a bismark.cov DNA methylation files into .fasta format by binerizing DNA methylation

################ How to run the program ################
#1. Use it as bash command
#2. Each sample should have both a bismark.cov file generated sing bismark_methylation_extraction
#3. The targets file needs to have the first 3 columns as 1) sample name, 2) tissue and 3) type

# The sript can be run as follows
# bismark_to_fasta.sh $odir $targets $cov_threshold $genome $genome_bismark $bed_location $bin_threshold $GATK_JAR $PICARD_JAR
########################################################

################## 
## Housekeeping ##
##################
# Check if necessary software is installed
type parallel >/dev/null 2>&1 || { echo >&2 "parallel is required but not installed.  Aborting."; exit 1; }
type samtools >/dev/null 2>&1 || { echo >&2 "samtools is required but not installed.  Aborting."; exit 1; }
type $GATK_JAR >/dev/null 2>&1 || { echo >&2 "gatk is required but not installed.  Aborting."; exit 1; }
type bedtools >/dev/null 2>&1 || { echo >&2 "bedtools is required but not installed.  Aborting."; exit 1; }
hash $PICARD_JAR >/dev/null 2>&1 || { echo >&2 "picard_JAR is required but not installed.  Aborting."; exit 1; }

# Check that all values have been supplied
if [[ $# -eq 0 ]] ; then
    echo 'No arguments supplied. Aborting.'
    exit 1
fi
if [ "$#" -ne 9 ]; then
  echo "The 9 required arguments not supplied. Aborting." >&2
  exit 1
fi

# Assign variables
if [ -z "$PS1" ]; then
    odir=$1 # ouput directory
    targets=$2 # phenotype dataframe
    cov_threshold=$3 # threshold for coverage to include cytosine in analysis
    genome=$4 # genome .fa file
    genome_bismark=$5
    bed_location=$6
    bin_threshold=$7
    GATK_JAR=$8
    PICARD_JAR=$9
fi

export njobs=4 #set this to a number so that minerva is happy

echo ""
echo "[ INPUTS ]"
echo ""
echo "output directory == $1"
echo "targets file == $2"
echo "coverage threshold == $3"
echo "FASTA file == $4"
echo "Location of bismark prepared .fa == $5"
echo "Location of bismark extracted methylation bismark.cov files == $6"
echo "threshold for binarization == $7"

# ensure targets file properly formatted
perl -pi -e 's/\r\n/\n/g' $targets

#################### run script ########################
# Apply coverage threshold to reduce bismark.cov file and convert to bed file
echo ""
echo "[ Checkiing all bismark.cov files are present ]"
echo ""

# ensure .bed files are present and make lists for sample types "interest" & "contrast"
rm bed_list
for file in $(awk '$3 ~ /interest|contrast/ {print $1}' $targets); do
    i=$bed_location/${file}.R1.fastq_paired_bismark_bt2_pe.bismark.cov.gz
    ls -d ${i} >> bed_list
if [ ! -f $(ls -d $bed_location/${file}.R1.fastq_paired_bismark_bt2_pe.bismark.cov.gz) ]
then
    echo "$0: File '${i}' not found." >&2
  exit 1
fi
done

# Apply coverage threshold to reduce bismark.cov file and convert to bed file
echo ""
echo "[ Applying coverage threshold to reduce reduce bismark.cov file and convert to sorted bedGraph file ]"
echo ""

function bismarkcov_to_reduced_bedGraph() {
    echo $1
    zcat $1 | awk '{ print $1"\t"($2-1)"\t"$3"\t"$4"\t"$5+$6;}' | awk -v a="$2" '$5 >= a {print $1"\t"$2"\t"$3"\t"$4}' | sort -k 1,1 -k2,2n | gzip > ${1%%.cov.gz}_reduced.bedGraph.gz
    echo ${1%%.cov.gz}_reduced.bedGraph.gz
}
export -f bismarkcov_to_reduced_bedGraph
parallel --xapply --will-cite -j+0 --eta bismarkcov_to_reduced_bedGraph {1} $cov_threshold :::: $odir/bed_list

# make output directory for each tissue and create list to pass to downstream functions
echo ""
echo "[ Making output directory for each tissue and create list to pass to downstream functions ]"
echo ""

rm $odir/input_bu
for file in $(awk '$3 ~ /interest|contrast/ {print $2}' $targets | uniq); do
    mkdir $odir/output_$file
    echo $odir/output_$file >> $odir/input_bu
done

# create symlinks to each samples combined_ms and combined_ums bedGraph files
echo ""
echo "[ Creating symlinks to each tissues BED files ]"
echo ""

grep 'interest\|contrast' $targets | while read line ; do
    set $line
    i=$bed_location/${1}.R1.fastq_paired_bismark_bt2_pe.bismark_reduced.bedGraph.gz
    echo $i
    ln -s $i $odir/output_${2}/${1}.R1.fastq_paired_bismark_bt2_pe.bismark_reduced.bedGraph.gz
done 

# merge bedgraphs from multiple samples into tissue bedGraph using bedtools unionbedg
echo ""
echo "[ Merging bedgraphs from multiple samples into tissue bedGraph using bedtools unionbedg ]"
echo ""

parallel --xapply --will-cite -j+0 --eta 'bedtools unionbedg -i {1}/*.R1.fastq_paired_bismark_bt2_pe.bismark_reduced.bedGraph.gz -filler NA > {1}/tissue.bedGraph' :::: $odir/input_bu

# sort tissue bedGraph using bedtools sort
echo ""
echo "[ Sorting tissue bedGraph using bedtools sort]"
echo ""

#parallel --xapply --will-cite -j+0 --eta 'bedtools sort -i {1}/tissue.bedGraph > {1}/tissue_sorted.bedGraph' :::: $odir/input_bu
parallel --xapply --will-cite -j+0 --eta 'sort -k 1,1 -k2,2n {1}/tissue.bedGraph > {1}/tissue_sorted.bedGraph' :::: $odir/input_bu

# Calculating the mean  DNA methylation accross all samples of tissue
echo ""
echo "[ Calculating mean methylation accross all samples of tissue ]"
echo ""

function mean_meth() { 
    awk '{sum=cnt=0; for (i=4;i<=NF;i++) if ($i != "NA") { sum+=$i; cnt++ } print $1"\t"$2"\t"$3"\t"(cnt ? sum/cnt : "NA") }' $1/tissue_sorted.bedGraph > $1/tissue_mc.bedGraph 
}
export -f mean_meth
parallel --xapply --will-cite -j+0 --eta mean_meth {1} :::: $odir/input_bu

# get strand info using bedtools
echo ""
echo "[ Obtaining strand information using bedtools ]"
echo ""

parallel --xapply -j+0 --eta 'bedtools getfasta -fi '$genome' -bed {1}/tissue_mc.bedGraph -tab -fo {1}/strand_info' :::: $odir/input_bu

# split to forward and reverse and label chr as in CT and GA bisulfite converted genomes
echo ""
echo "[ Splitting to forward and reverse and label chr as in CT and GA bisulfite converted genomes ]"
echo ""

function split_bed() {
    bin_threshold=$2
paste $1/tissue_mc.bedGraph $1/strand_info | awk '$6=="C";$6=="c" {print $1, $2, $3, $4, $5}' OFS='\t' > $1/temp2.bed
awk -v a="${bin_threshold}" '{print $1"_CT_converted",$3,$3,($4>=a)? "T" : ($4<a)? "C" : $4, ($4>=a)? "C" : ($4<a)? "T" : $4}' OFS='\t' $1/temp2.bed | sort -k1,1 -k2,2n > $1/C2T.bed
paste $1/tissue_mc.bedGraph $1/strand_info | awk '$6=="G";$6=="g" {print $1, $2, $3, $4, $5}' OFS='\t' > $1/temp3.bed
awk -v a="${bin_threshold}" '{print $1"_GA_converted",$3,$3,($4>=a)? "A" : ($4<a)? "G" : $4, ($4>=a)? "G" : ($4<a)? "A" : $4}' OFS='\t' $1/temp3.bed | sort -k1,1 -k2,2n > $1/G2A.bed
}
export -f split_bed

parallel --xapply --will-cite -j+0 --eta split_bed {1} $bin_threshold :::: $odir/input_bu

### convert to .vcf ###
# create .vcf header
echo ""
echo "[ Creating .vcf header ]"
echo ""

echo "##fileformat=VCFv4.1" > $odir/vcf_header
echo "##FILTER=<ID=LowQual,Description=" >> $odir/vcf_header
echo "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=" >> $odir/vcf_header
echo "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=" >> $odir/vcf_header
echo "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=" >> $odir/vcf_header
echo "##FORMAT=<ID=GT,Number=1,Type=String,Description=" >> $odir/vcf_header
echo "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=" >> $odir/vcf_header
echo "##GATKCommandLine.HaplotypeCaller=<ID=HaplotypeCaller,Version=3.4-3-gd1ac142,Date=" >> $odir/vcf_header
echo "##INFO=<ID=AC,Number=A,Type=Integer,Description=" >> $odir/vcf_header
echo "##INFO=<ID=AF,Number=A,Type=Float,Description=" >> $odir/vcf_header
echo "##INFO=<ID=AN,Number=1,Type=Integer,Description=" >> $odir/vcf_header
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTISSUE" >> $odir/vcf_header
echo -e "100\tPASS\t.\tGT:AD:DP:GQ:PL\t1/1:141:282:99:255,0,25" > $odir/redundent_vcf_vals # add redundent info to aditional fields of .vcf

function bed_to_vcf() {
    NUM_LINES=$(wc -l $1/C2T.bed | awk '{print $1}')
    paste $1/C2T.bed <(awk -v a="$NUM_LINES" '{for(i=0;i<a;i++)print}' redundent_vcf_vals) > $1/temp1.vcf
    cat ./vcf_header $1/temp1.vcf > $1/methylotypes_CT.vcf
    rm $1/temp1.vcf

    NUM_LINES=$(wc -l $1/G2A.bed | awk '{print $1}')
    paste $1/G2A.bed <(awk -v a="$NUM_LINES" '{for(i=0;i<a;i++)print}' redundent_vcf_vals) > $1/temp2.vcf
    cat ./vcf_header $1/temp2.vcf > $1/methylotypes_GA.vcf
    rm  $1/temp2.vcf
}
export -f bed_to_vcf

parallel --xapply --will-cite -j+0 --eta bed_to_vcf {1} :::: $odir/input_bu 
rm ./vcf_header ./redundent_vcf_vals

# define function for vcf2fasta 
function vcf2fasta_wholegenome() {
# directory of bisulfite genome prepared using genome_prepare.sh (implements bismark genome preparation with dictionary)
FADIR=$1
# directory of .vcf file & output directory
VCF_DIR=$2
OUTDIR=$2

CT_vcf=$2/methylotypes_CT.vcf
GA_vcf=$2/methylotypes_GA.vcf

PICARD_JAR=$3
GATK_JAR=$4

## CT_genome ##
# sort methylotype.vcf in same order as .fasta
rm $OUTDIR/methylotypes_CT.sorted.vcf
$PICARD_JAR SortVcf I=$CT_vcf O=$OUTDIR/methylotypes_CT.sorted.vcf SEQUENCE_DICTIONARY=$FADIR/CT_conversion/genome_mfa.CT_conversion.dict 

# create methylotype specific .fasta files
rm $OUTDIR/methylotypes_CT.sorted.vcf.idx
line_width=50
$GATK_JAR -T FastaAlternateReferenceMaker -R $FADIR/CT_conversion/genome_mfa.CT_conversion.fa -o $OUTDIR/methylotype.CT_conversion.fasta -V $OUTDIR/methylotypes_CT.sorted.vcf -lw $line_width

## GA_genome ##
# sort methylotype.vcf in same order as .fasta
rm $OUTDIR/methylotypes_GA.sorted.vcf
$PICARD_JAR SortVcf I=$GA_vcf O=$OUTDIR/methylotypes_GA.sorted.vcf SEQUENCE_DICTIONARY=$FADIR/GA_conversion/genome_mfa.GA_conversion.dict 

# # create methylotype specific .fasta files
rm $OUTDIR/methylotypes_GA.sorted.vcf.idx
line_width=50
$GATK_JAR -T FastaAlternateReferenceMaker -R $FADIR/GA_conversion/genome_mfa.GA_conversion.fa -o $OUTDIR/methylotype.GA_conversion.fasta -V $OUTDIR/methylotypes_GA.sorted.vcf -lw $line_width

# format and combine CT and GA fasta files
sed -e 's/\(>\).*\(\ \)/\1\2/' $OUTDIR/methylotype.CT_conversion.fasta | tr -d "[:blank:]" | sed -e 's/converted:1/converted/g' > $OUTDIR/temp.fasta
mv $OUTDIR/temp.fasta $OUTDIR/methylotype.CT_conversion.fasta

sed -e 's/\(>\).*\(\ \)/\1\2/' $OUTDIR/methylotype.GA_conversion.fasta | tr -d "[:blank:]" | sed -e 's/converted:1/converted/g' > $OUTDIR/temp.fasta
mv $OUTDIR/temp.fasta $OUTDIR/methylotype.GA_conversion.fasta
}
export -f vcf2fasta_wholegenome

# run vcf2fasta in parallel
echo ""
echo "[ Converting VCF to FASTA ]"
echo ""
parallel --xapply --will-cite -j+0 --eta vcf2fasta_wholegenome $genome_bismark {1} $PICARD_JAR $GATK_JAR :::: $odir/input_bu

### break-up fasta files based on coverage of all tissues ###
# get combined coverage .bed file
rm $odir/intersect_C2T.bed $odir/intersect_G2A.bed
multiIntersectBed -i output_*/C2T.bed >> $odir/intersect_C2T.bed
multiIntersectBed -i output_*/G2A.bed >> $odir/intersect_G2A.bed

# sort combined .bed files
rm $odir/intersect_C2T_sorted.bed $odir/intersect_G2A_sorted.bed
sort -k 1,1 -k2,2n $odir/intersect_C2T.bed >> $odir/intersect_C2T_sorted.bed 
sort -k 1,1 -k2,2n $odir/intersect_G2A.bed >> $odir/intersect_G2A_sorted.bed 

# remove lines without coverage in all bedGraphs
ntiss=$(find output_* -type d | wc -l)
awk -v a="$ntiss" '$4 != a {print $1"\t"$2"\t"$3"\t"$4}' $odir/intersect_C2T_sorted.bed > $odir/intersect_C2T_sub.bed
awk -v a="$ntiss" '$4 != a {print $1"\t"$2"\t"$3"\t"$4}' $odir/intersect_G2A_sorted.bed > $odir/intersect_G2A_sub.bed
#rm $odir/intersect_C2T.bed $odir/intersect_G2A.bed

# create intersect with C2T/G2A genome.bed files
if [ -s $odir/intersect_C2T_sub.bed ]; then
    bedtools complement -i $odir/intersect_C2T_sub.bed -g $genome_bismark/CT_conversion/genome_mfa.CT_conversion.genome  > $odir/genome_intersect_CT.bed
else
    awk '{print $1"\t"'0'"\t"$2}' $genome_bismark/CT_conversion/genome_mfa.CT_conversion.genome  > $odir/genome_intersect_CT.bed
fi

if [ -s $odir/intersect_G2A_sub.bed ]; then
    bedtools complement -i $odir/intersect_G2A_sub.bed -g $genome_bismark/GA_conversion/genome_mfa.GA_conversion.genome  > $odir/genome_intersect_GA.bed
else
    awk '{print $1"\t"'0'"\t"$2}' $genome_bismark/GA_conversion/genome_mfa.GA_conversion.genome  > $odir/genome_intersect_GA.bed
fi

# break-up .fasta files 
parallel --xapply -j+0 --eta 'bedtools getfasta -fi {1}/methylotype.CT_conversion.fasta -bed '$odir/genome_intersect_CT.bed' -fo {1}/temp.fasta' :::: $odir/input_bu
parallel --xapply -j+0 --eta mv {1}/temp.fasta {1}/methylotype.CT_conversion.fasta :::: $odir/input_bu

parallel --xapply -j+0 --eta 'bedtools getfasta -fi {1}/methylotype.GA_conversion.fasta -bed '$odir/genome_intersect_GA.bed' -fo {1}/temp.fasta' :::: $odir/input_bu
parallel --xapply -j+0 --eta mv {1}/temp.fasta {1}/methylotype.GA_conversion.fasta :::: $odir/input_bu

### rename "chromosomes" to include tissue name and combine all into master for indexing ###
for i in $(cat $odir/input_bu); do
    echo $i
    temp=${i#*output_}
    echo $temp
    sed -i "/>/ s/$/:methylotype_${temp%/*}/g" ${i}/methylotype.CT_conversion.fasta
    sed -i "/>/ s/$/:methylotype_${temp%/*}/g" ${i}/methylotype.GA_conversion.fasta
done

cat $odir/output_*/*conversion.fasta > master_methylotype.fasta

########### end script ###########
echo "[ END - bedGraph_to_fasta ]"
