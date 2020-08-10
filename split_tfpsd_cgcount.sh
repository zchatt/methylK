#!/bin/bash
# split true and false uniquely pseudoaligned reads and count the cg content

###########
## INPUT ##
###########

SAMPLE_NAME=$1
OUTDIR=$2
mkdir $OUTDIR
SAM=$3
TISSUE=$4

#########
## RUN ##
#########

# define cg count function
function count_cg() {
    ASSAY=$1
    OUTDIR=$2

rm $OUTDIR/${ASSAY}_counts2 $OUTDIR/${ASSAY}_c_occurance $OUTDIR/${ASSAY}_g_occurance $OUTDIR/${ASSAY}_tempC* $OUTDIR/${ASSAY}_tempG* 2> /dev/null
# subset based on assay
grep $ASSAY $OUTDIR/temp_counts > $OUTDIR/${ASSAY}_counts2
# count occurances of C ($2)
awk '{print $2}' $OUTDIR/${ASSAY}_counts2 | sort | uniq -c | awk '{ print $2 " " $1}' | sort -n > $OUTDIR/${ASSAY}_c_occurance
# count occurances of G ($3)
awk '{print $3}' $OUTDIR/${ASSAY}_counts2 | sort | uniq -c | awk '{ print $2 " " $1}' | sort -n > $OUTDIR/${ASSAY}_g_occurance
# joing the C and G counts for the length of the read and assign to assay
join -e- -a1 -a2 $OUTDIR/${ASSAY}_c_occurance -o 0 1.2 $OUTDIR/temp.seq | column -t > $OUTDIR/${ASSAY}_tempC2
awk '!x[$1]++' $OUTDIR/${ASSAY}_tempC2 | sort -n > $OUTDIR/${ASSAY}_tempC
join -e- -a1 -a2 $OUTDIR/${ASSAY}_g_occurance -o 0 1.2 $OUTDIR/temp.seq | column -t > $OUTDIR/${ASSAY}_tempG2
awk '!x[$1]++' $OUTDIR/${ASSAY}_tempG2 | sort -n > $OUTDIR/${ASSAY}_tempG
join -e- -a1 -a2 $OUTDIR/${ASSAY}_tempC -o 0 1.2 $OUTDIR/${ASSAY}_tempG -o 2.2 | column -t | sort -n > $OUTDIR/${ASSAY}_cgcounttemp

# calculate the proportion of each co-methylated event of all pseudoaligned reads to assay, regardless of tissue assignment
rm $OUTDIR/${ASSAY}_counts2 $OUTDIR/${ASSAY}_c_occurance $OUTDIR/${ASSAY}_g_occurance $OUTDIR/${ASSAY}_tempC* $OUTDIR/${ASSAY}_tempG* 2> /dev/null
}
export -f count_cg


# select uniquely pseudoaligned reads to one cells methylotpye from .sam file output
awk '$12 == "NH:i:1" {print $0}' $SAM | sort > $OUTDIR/${SAMPLE_NAME}_upsd.sam

# combine the sequences of paired reads fo cg counting
awk '{if(a!=$1) {a=$1; printf "\n%s%s",$10,FS} else {a=$1;$1="";printf $10 }} END {printf "\n" }' $OUTDIR/${SAMPLE_NAME}_upsd.sam > $OUTDIR/tmp1
awk '{if(a!=$1) {a=$1; printf "\n%s%s",$0,FS} else {a=$1;$1="";printf $1 }} END {printf "\n" }' $OUTDIR/${SAMPLE_NAME}_upsd.sam > $OUTDIR/tmp2
paste -d'\0' <(awk '{print $1}' $OUTDIR/tmp1) <(awk '{print $2}' $OUTDIR/tmp1) > $OUTDIR/tmp4
paste <(awk '{print $1" "$3}' $OUTDIR/tmp2) $OUTDIR/tmp4 -d "  " > $OUTDIR/tmp5

# select true and false
grep $TISSUE $OUTDIR/tmp5 | awk 'NF > 0' > $OUTDIR/${SAMPLE_NAME}_true_psd.sam
grep -v $TISSUE $OUTDIR/tmp5 | awk 'NF > 0' > $OUTDIR/${SAMPLE_NAME}_false_psd.sam

# determine read length from first 1000 reads
chars=$(head -1000 $OUTDIR/${SAMPLE_NAME}_true_psd.sam | awk '{print $3}' | wc -c)
words=$(head -1000 $OUTDIR/${SAMPLE_NAME}_true_psd.sam | awk '{print $3}' | wc -w)
read_length=$(( ${chars} / ${words} ))

# makes sequence of read length
seq 0 $read_length | sort -n > $OUTDIR/temp.seq

#####################################################
### counts for true (matched pseudoaligend reads) ###
#####################################################
SAM_M=$OUTDIR/${SAMPLE_NAME}_true_psd.sam

# count occurances of "C" within each "CT" assay and "G" within each "GA" assay
awk '{print $3}' $SAM_M | awk '{print $2, gsub("C","")}' > $OUTDIR/c_counts
awk '{print $3}' $SAM_M | awk '{print $2, gsub("G","")}' > $OUTDIR/g_counts
awk '{print $2}' $SAM_M > $OUTDIR/temp_assay
paste $OUTDIR/temp_assay $OUTDIR/c_counts $OUTDIR/g_counts > $OUTDIR/temp_counts2
sort $OUTDIR/temp_counts2 > $OUTDIR/temp_counts

awk '!seen[$1]++' $OUTDIR/temp_counts | awk '{print $1}' | sort > $OUTDIR/assay_list
parallel --xapply --will-cite -j+0 --eta count_cg {1} $OUTDIR :::: $OUTDIR/assay_list

# combine all counts into single file of sample
for filename2 in $(ls $OUTDIR/*_cgcounttemp)
do
filename=${filename2##*/}
echo -e "counts ${filename//cgcounttemp/C} ${filename//cgcounttemp/G}" | cat - $filename2 > ${filename2}.new
#echo -e "counts ${filename//cgcounttemp/C} ${filename//cgcounttemp/G}"
#echo Done ${filename2}
done
#  join all files and remove temp files
paste $OUTDIR/*cgcounttemp.new > $OUTDIR/${SAMPLE_NAME}_true_cgcount

# remove temp cout files
rm $OUTDIR/*_cgcounttemp $OUTDIR/*_cgcounttemp.new $OUTDIR/c_counts $OUTDIR/g_counts $OUTDIR/temp_assay $OUTDIR/temp_counts $OUTDIR/temp_counts2 $OUTDIR/tmp{1,2,3,4,5} $OUTDIR/assay_list

##########################################################
### counts for false (non-matched pseudoaligend reads) ###
##########################################################
SAM_UM=$OUTDIR/${SAMPLE_NAME}_false_psd.sam

# count occurances of "C" within each "CT" assay and "G" within each "GA" assay
awk '{print $3}' $SAM_UM | awk '{print $2, gsub("C","")}' > $OUTDIR/c_counts
awk '{print $3}' $SAM_UM | awk '{print $2, gsub("G","")}' > $OUTDIR/g_counts
awk '{print $2}' $SAM_UM > $OUTDIR/temp_assay
paste $OUTDIR/temp_assay $OUTDIR/c_counts $OUTDIR/g_counts > $OUTDIR/temp_counts2
sort $OUTDIR/temp_counts2 > $OUTDIR/temp_counts

awk '!seen[$1]++' $OUTDIR/temp_counts | awk '{print $1}' | sort > $OUTDIR/assay_list
parallel --xapply --will-cite -j+0 --eta count_cg {1} $OUTDIR :::: $OUTDIR/assay_list 2> /dev/null

# combine all counts into single file of sample
for filename2 in $(ls $OUTDIR/*_cgcounttemp)
do
filename=${filename2##*/}
echo -e "counts ${filename//cgcounttemp/C} ${filename//cgcounttemp/G}" | cat - $filename2 > ${filename2}.new
#echo -e "counts ${filename//cgcounttemp/C} ${filename//cgcounttemp/G}"
#echo Done ${filename2}
done
#  join all files and remove temp files
paste $OUTDIR/*cgcounttemp.new > $OUTDIR/${SAMPLE_NAME}_false_cgcount

# remove temporary files
rm $OUTDIR/c_counts $OUTDIR/g_counts $OUTDIR/temp_assay $OUTDIR/temp_counts $OUTDIR/temp_counts2 $OUTDIR/*_proptemp $OUTDIR/*proptemp.new $OUTDIR/assay_list 2> /dev/null
rm $OUTDIR/tmp.sam $OUTDIR/temp* $OUTDIR/c_* $OUTDIR/g_* $OUTDIR/*_cgcounttemp $OUTDIR/*cgcounttemp.new $OUTDIR/tempC* $OUTDIR/tempG* 2> /dev/null
rm $OUTDIR/sam* $OUTDIR/*_tmp $OUTDIR/*_tmp3 $OUTDIR/*_tmp2 $OUTDIR/tmp* 2> /dev/null