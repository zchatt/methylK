#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
stringsAsFactors=FALSE

# check if all arguments present
if (length(args)==0) {
  stop("At least one argument must be supplied", call.=FALSE)
} else if (length(args) < 3) {
  # default output file
  stop("The 3 input variables have not been supplied", call.=FALSE)
}

# Converting possible relative path to absolute path for the other arguments
for (i in c(2)){
  args[i]<-scan(pipe(paste("readlink -m ", args[i])),"",quiet=TRUE)
}

### ZC testing ###
# methylK_dir="/project/RDS-SMS-FFbigdata-RW/Epigenetics/tngbs/zacchatterton_tNGBS060219/mk_organise"
# working_dir="/project/RDS-SMS-FFbigdata-RW/Epigenetics/tngbs/zacchatterton_tNGBS060219/mk_organise/test/output"
# targets<-read.delim(file="/project/RDS-SMS-FFbigdata-RW/Epigenetics/tngbs/zacchatterton_tNGBS060219/mk_organise/test/targets.txt",sep='\t', header=T) 

# module load R/3.6.3
# Rscript --vanilla $methylK_dir/snr_calc.R $odir $methylK_dir $targets
###         ###

#############
### INPUT ###
#############

working_dir=args[1]
methylK_dir=args[2]
snr_thresh_path=args[3]
targets=read.delim(file=args[4],sep='\t', header=T) #targets file should have column with headers 1) "sample_name", 2) "tissue" , 3) "counts_fq"
source(paste0(methylK_dir,"/","read_cgcnt.R"))

#############################
### LIBRARY AND FUNCTIONS ###
#############################
library(plyr)
library(dplyr)
library(tidyr)

############################
### READ DATA AND FORMAT ###
############################
# read data
setwd(working_dir)
cg_count=read_cgcnt(working_dir)
tiss_interest<-unique(as.character(targets$tissue[targets$type == "interest"]))
targets<-targets[targets$type %in% c("identify"),]

# calculate the number of reads within the trimmed, paired and truncated .fastq files for each sample types "interest" or "contrast"
for (i in 1:nrow(targets)){
  sample_n=targets$sample_seqname[i]
  res=system(paste0("zcat ",sample_n,".paired_truncated_R2.fastq.gz | wc -l", sep=""),intern=TRUE)
  targets$counts_fq[i]<-as.numeric(unlist(strsplit(res," "))[1])/4
}

##############
### FORMAT ###
##############
# get assay details
df<-cg_count
df$assay<-gsub(".methylotype_.*._","",df$assay_long)

# map tissue and number of fq reads
sample_tissue<- mapvalues(as.character(df$sample_name), from=as.character(targets$sample_seqname), to=as.character(targets$tissue))
counts_fq<- mapvalues(as.character(df$sample_name), from=as.character(targets$sample_seqname), to=as.character(targets$counts_fq))
df<-cbind(df,sample_tissue,counts_fq)

# calculate fractional values of unique pseudo / total reads for each sample
df$frac_upsdfq<-as.numeric(as.character(df$u_psdcnt))/as.numeric(as.character(df$counts_fq))

###########################
### THRESHOLD USING SNR ###
###########################

df = df %>% complete(sample_name, nesting(assay_long,counts,assay,psd_cell))
col_interest<-c("counts","assay_long","sample_name","assay","psd_cell","frac_upsdfq","frac_upsdfq_thresh")
for (z in 1:length(tiss_interest)){
  tryCatch({ 
    res<-c("counts","assay_long","sample_name","assay","psd_cell","frac_upsdfq","frac_upsdfq_thresh")
    cell1<-tiss_interest[z]
    
    # COI file name
    file_name=list.files(path = snr_thresh_path)[grep(paste0(tiss_interest[z],sep="_thresh-"),list.files(path = snr_thresh_path))]
    
    # read in SNR threshold file
    snr_thresh<-read.table(file=paste0(snr_thresh_path,sep="/",file_name), sep='\t', header=T)
    
    # subset df using threshold assays
    mvals<-paste(df$assay_long,df$counts)
    
    aval<-paste(snr_thresh$assay_long,snr_thresh$cometh_count)
    df_thresh<-df[mvals %in% aval,]
    
    # fill in missing data/ sample
    d1= df_thresh %>% complete(counts, nesting(assay_long,sample_name,assay,psd_cell))
    
    # subtract SNR from each assay_long
    snr<-1/as.numeric(mapvalues(as.character(d1$assay_long), from=as.character(snr_thresh$assay_long), to=as.character(snr_thresh$SNR)))
    d1$frac_upsdfq_thresh<-d1$frac_upsdfq-snr
    d1$frac_upsdfq_thresh[d1$frac_upsdfq_thresh < 0]<-0
    
    # subset and bind to res
    d1<-d1[,col_interest]
    res<-rbind(res,d1)
    
    # write results to file
    res<-res[-1,]
    write.table(res,file=paste0(paste(tiss_interest[z], sep="", collapse=""),"_mkdf.txt"), sep='\t', quote = FALSE)

  },
  error=function(e){})
}
