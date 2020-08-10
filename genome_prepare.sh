#!/bin/bash
# The following scripts create a bisulfite genome from .fasta, index and create a dictionary necessary for .vcf to .fasta conversion
# reference genome directory
BISGENOME=$1

####################################
### Bisulfite Genome Preparation ###
####################################
# bismark
bismark_genome_preparation $BISGENOME
# bwa-meth
#bwameth.py index $BISGENOME/*.fa

#######################################################
### Index Bisulfite Genome & Create .bed for Genome ###
#######################################################

FADIR=$BISGENOME/Bisulfite_Genome
samtools faidx $FADIR/CT_conversion/genome_mfa.CT_conversion.fa
samtools faidx $FADIR/GA_conversion/genome_mfa.GA_conversion.fa

cut -f1,2 $FADIR/CT_conversion/genome_mfa.CT_conversion.fa.fai > $FADIR/CT_conversion/genome_mfa.CT_conversion.genome
cut -f1,2 $FADIR/GA_conversion/genome_mfa.GA_conversion.fa.fai > $FADIR/GA_conversion/genome_mfa.GA_conversion.genome

##########################################
### Create Bisulfite Genome Dictionary ###
##########################################
function picard_else
{
    if type picard >/dev/null
    then  
    picard CreateSequenceDictionary R=$FADIR/CT_conversion/genome_mfa.CT_conversion.fa O=$FADIR/CT_conversion/genome_mfa.CT_conversion.dict
	picard CreateSequenceDictionary R=$FADIR/GA_conversion/genome_mfa.GA_conversion.fa O=$FADIR/GA_conversion/genome_mfa.GA_conversion.dict

    else  
    java -jar $PICARD_JAR CreateSequenceDictionary R=$FADIR/CT_conversion/genome_mfa.CT_conversion.fa O=$FADIR/CT_conversion/genome_mfa.CT_conversion.dict
	java -jar $PICARD_JAR CreateSequenceDictionary R=$FADIR/GA_conversion/genome_mfa.GA_conversion.fa O=$FADIR/GA_conversion/genome_mfa.GA_conversion.dict
    fi
}
picard_else

