#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
stringsAsFactors=FALSE

#############
### INPUT ###
#############
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
working_dir=args[1]
methylK_dir=args[2]
targets=read.delim(file=args[3],sep='\t', header=T) #targets file should have column with headers 1) "sample_name", 2) "tissue" , 3) "counts_fq"

targets<-targets[targets$type %in% c("interest","contrast"),]
tiss_interest<-unique(as.character(targets$tissue[targets$type == "interest"]))
tiss_contrast<-unique(as.character(targets$tissue[targets$type == "contrast"]))


### ZC testing ###
#methylK_dir="/project/RDS-SMS-FFbigdata-RW/Epigenetics/tngbs/zacchatterton_tNGBS060219/mk_organise"
#working_dir="/project/RDS-SMS-FFbigdata-RW/Epigenetics/tngbs/zacchatterton_tNGBS060219/mk_organise/test/output"
#df_psd<-format_psd_cell("/project/RDS-SMS-FFbigdata-RW/Epigenetics/tngbs/zacchatterton_tNGBS060219/mk_organise/test/output")
#targets<-read.delim(file="/project/RDS-SMS-FFbigdata-RW/Epigenetics/tngbs/zacchatterton_tNGBS060219/mk_organise/test/targets.txt",sep='\t', header=T) 
#targets<-targets[targets$type %in% c("interest","contrast"),]

# module load R/3.6.3
# Rscript --vanilla $methylK_dir/snr_calc.R $odir $methylK_dir $targets
###         ###


#############################
### LIBRARY AND FUNCTIONS ###
#############################
source(paste0(methylK_dir,"/","read_cgcnt.R"))
source(paste0(methylK_dir,"/","format-psdcnt.R"))
library(plyr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(scales)
library(caret)

############################
### READ DATA AND FORMAT ###
############################
# read data
setwd(working_dir)
cg_count=read_cgcnt(working_dir)
df_psd<-format_psd_cell(working_dir)

# calculate the number of reads within the trimmed, paired and truncated .fastq files for each sample types "interest" or "contrast"
for (i in 1:nrow(targets)){
  sample_n=targets$sample_seqname[i]
res=system(paste0("zcat ",sample_n,".paired_truncated_R2.fastq.gz | wc -l", sep=""),intern=TRUE)
targets$counts_fq[i]<-as.numeric(unlist(strsplit(res," "))[1])/4
}

# calculate true/false ratio for each assay/cgcount
tf<-1:nrow(cg_count) %in% grep("_true_cgcount",cg_count$sample)
df<-cbind(cg_count,tf)
df$assay<-gsub(".methylotype_.*._","",df$assay_long)
df$sample_name<-gsub("_false","",df$sample_name)
df$sample_name<-gsub("_true","",df$sample_name)

# map tissue and number of fq reads
sample_tissue<- mapvalues(as.character(df$sample_name), from=as.character(targets$sample_seqname), to=as.character(targets$tissue))
counts_fq<- mapvalues(as.character(df$sample_name), from=as.character(targets$sample_seqname), to=as.character(targets$counts_fq))
df<-cbind(df,sample_tissue,counts_fq)

# calculate fractional values of unique pseudo / total reads for each sample
df$frac_upsdfq<-as.numeric(as.character(df$u_psdcnt))/as.numeric(as.character(df$counts_fq))

# calculate fractional values of unique pseudo / total unique pseudo for each assay
df$frac_psdpsd<-df$counts
p<-unique(df$sample_name)
# for each sample
for (i in 1:length(unique(p))){
  a2<-df[df$sample_name == as.character(p[i]),]
# for each assay 
  ua<-unique(a2$assay_long)
  for (z in 1:length(unique(ua))){
    c2<-df$u_psdcnt[df$assay_long == as.character(ua[z]) & df$sample_name == as.character(p[i])]/sum(a2$u_psdcnt[a2$assay_long == ua[z] ],na.rm = T)
    df$frac_psdpsd[df$assay_long == as.character(ua[z]) & df$sample_name == as.character(p[i])]<-c2
  }}

#######################
### Calculating SNR ###
#######################

# make header to bind 
snr_total<-c("assay_long","cometh_count","TPR","FPR","SNR")

# map total counts to each assay
d1a<-aggregate(count ~ assay * sample_name,df_psd,function(x) sum(x,na.rm = TRUE))
map_val1<-paste0(d1a$sample_name,d1a$assay)
d1a<-cbind(d1a,map_val1)

a1<-gsub("_CT.*","",df$assay)
a1<-gsub("_GA.*","",a1)
map_val2<-paste0(df$sample_name,a1)
counts_total_assay<- mapvalues(as.character(map_val2), from=as.character(d1a$map_val1), to=as.character(d1a$count))
df<-cbind(df,counts_total_assay)

# loop through COI's
gl<-list()
for (i in 1:length(tiss_interest)){
  cell1<-as.character(tiss_interest[i])
  
  # calculate tpr and for for COI
  d2a<-df[df$sample_tissue == cell1 & df$tf == "TRUE",]
  tpr1<-as.numeric(as.character(d2a$u_psdcnt))/
    as.numeric(as.character(d2a$counts_total_assay))
  
  d2b<-df[df$sample_tissue %in% tiss_contrast & df$tf == "FALSE",]
  fpr1<-as.numeric(as.character(d2b$u_psdcnt))/
    as.numeric(as.character(d2b$counts_total_assay))
  
  d2a<-cbind(d2a,tpr1)
  d2b<-cbind(d2b,fpr1)
  
  # aggregate  - calc mean TPR and FPR for each assay_long + co-meth + tissue
  TPR_ag<-aggregate(tpr1 ~ assay_long * counts, d2a,function(x) mean(x,na.rm = TRUE))
  FPR_ag<-aggregate(fpr1 ~ assay_long * counts, d2b,function(x) mean(x,na.rm = TRUE))
  
  snr_pre<-merge(TPR_ag,FPR_ag,by=c("assay_long","counts"), all.x = TRUE)
  
  snr_pre$tpr1[is.na(snr_pre$tpr1)] <- 0
  snr_pre$fpr1[is.na(snr_pre$fpr1)] <- 0
  snr<-cbind(snr_pre,snr=snr_pre$tpr1/snr_pre$fpr1)
  snr$snr[snr_pre$fpr1 == 0] <- 1/(snr_pre$tpr1[snr_pre$fpr1 == 0]) # if pbmc = 0 we invert the coi signal to define the signal-to-noise
  
  colnames(snr)<-c("assay_long","cometh_count","TPR","FPR","SNR")
  write.table(snr,file=paste0(cell1,sep="_","snr.txt"), sep='\t')
  snr_total<-rbind(snr_total,snr)
  
}

#################################
### DEFINE SNR THRESHOLD/ COI ###
#################################

gl<-list()
# for each COI calculate SNR thresholds
for (z in 1:length(tiss_interest)){
  tryCatch({  
  # data frame & fill in missing with 0
  dfa<-df[df$sample_tissue %in% c(tiss_interest[z],tiss_contrast),]
  dfa$frac_upsdfq[is.na(dfa$frac_upsdfq)]<-0
  
  # read in SNR for COI
  snr_total<-read.table(file=paste0(tiss_interest[z],sep="_","snr.txt"), sep='\t', header=T)

res2<-0
thresholds<-seq(0,20000,10)
for (i2 in 1:length(thresholds)){
  # threshold
  thresh<-thresholds[i2]
  
  # cell-type
  #p<- as.character(unique(dfa$psd_cell))
  
  # get assays that pass snr
  a1<-snr_total[as.numeric(snr_total$SNR) >= thresh, ]
  aval<-paste(a1$assay_long,a1$cometh_count)
  
  # subset df using threshold assays
  mvals<-paste(dfa$assay_long,dfa$counts)
  dfa2<-dfa[mvals %in% aval,]
  
  # calc precision, recall and F1 for each cell-combination
  cells<-unique(c(tiss_interest[z],tiss_contrast))
  length(cells)
  comb_table<-t(combn(c(as.character(cells),rev(as.character(cells))),2))
  comb_table<-t(comb_table[!duplicated(comb_table),])
  res<-as.data.frame(matrix(0,ncol(comb_table),3))
  colnames(res)<-c("Precision","Recall","F1")
  for (i in 1:ncol(comb_table)) {
    
    cell1=comb_table[1,i]
    cell2=comb_table[2,i]
    # create contigency table
    tmp<-matrix(0,2,2)
    # number of COI == COI
    tmp[1,1]<-sum(as.numeric(dfa2$frac_upsdfq[dfa2$psd_cell == cell1 & dfa2$sample_tissue == cell1]))
    # number of COI == CONTRAST
    tmp[2,1]<-sum(as.numeric(dfa2$frac_upsdfq[dfa2$psd_cell == cell2 & dfa2$sample_tissue == cell1]))
    # number of CONTRAST == CONTRAST
    tmp[2,2]<-sum(as.numeric(dfa2$frac_upsdfq[dfa2$psd_cell == cell2 & dfa2$sample_tissue == cell2]))
    # number of CONTRAST == COI
    tmp[1,2]<-sum(as.numeric(dfa2$frac_upsdfq[dfa2$psd_cell == cell1 & dfa2$sample_tissue== cell2]))
    
    row.names(res)[i]<-paste(cell1,cell2,sep="-")
    res[i,1]<-precision(as.table(tmp))
    res[i,2]<-recall(as.table(tmp))
    res[i,3]<-F_meas(as.table(tmp))
    
  }
  res2<-c(res2,res[1,3])
}
F1<-res2[-1]

thresholds<-thresholds[!is.na(F1)] # added 160720
F1<-F1[!is.na(F1)]# added 160720

# calculate the inflection point i.e. where smoothed F1 changes sign
thrshold <- 0.0001
d1<-diff(F1)
inflection_pt<-which.max(d1 / max(d1) < thrshold)

lo <- loess(F1~thresholds)
xl <- seq(min(thresholds),max(thresholds), (max(thresholds) - min(thresholds))/1000)
out = predict(lo,xl)
out[out > 1]<-1
infl <- c(FALSE, diff(diff(out)>0)!=0)
inflection_pt<-which.max(infl)

# assign to plot
dplot<-as.data.frame(cbind(thresholds,F1))
gl[[z]]<-ggplot(dplot, aes(x=thresholds, y=F1)) + 
  geom_point()+ theme_minimal() + 
  ylab("F1") + xlab("SNR threshold") + ggtitle(tiss_interest[z]) + geom_vline(xintercept=thresholds[inflection_pt], col="red",lty=2) +
  geom_hline(yintercept=F1[inflection_pt], col="blue",lty=2) +
  labs(subtitle =paste0("SNR threshold =",sep=" ",thresholds[inflection_pt],"; F1 =",sep=" ",round(F1[inflection_pt],5)))

# threhold SNR file for COI
snr<-snr_total[as.numeric(snr_total$SNR) >= thresholds[inflection_pt], ]
write.table(snr,file=paste0(tiss_interest[z],sep="_thresh-",thresholds[inflection_pt],"-snr.txt"), sep='\t')
  },
error=function(e){})
}










