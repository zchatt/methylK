#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

####################
###  Run Example ###
####################

# WD="/Users/zacc/Documents/MtSinai/Codes/tNGBS/kmer_tngbs_codes/FTD_cfDNA/CX_reports/"
# assay="/Users/zacc/Documents/MtSinai/Codes/tNGBS/kmer_tngbs_codes/assay_description_small.txt"
# Rdata="/Users/zacc/Documents/MtSinai/Codes/tNGBS/kmer_tngbs_codes/FTD_cfDNA/CX_reports/tngbs_cfdna_dataframes.Rdata"

# Rscript --vanilla WFR_estimates.R $WD $assay $Rdata

#############
### INPUT ###
#############

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args) < 3) {
  # default output file
  stop("The 3 input variables hav not been supplied", call.=FALSE)
}
setwd(args[1])
assay=read.delim(file=args[2], sep='\t', header=TRUE)
load(args[3])

WD="/Users/zacc/Documents/MtSinai/Codes/tNGBS/kmer_tngbs_codes/FTD_cfDNA/CX_reports/"
assay="/Users/zacc/new_methylK/assay_description_small.txt"
Rdata="/Users/zacc/Documents/MtSinai/Codes/tNGBS/kmer_tngbs_codes/FTD_cfDNA/CX_reports/tngbs_cfdna_dataframes.Rdata"

setwd(WD)
assay=read.delim(file=assay, sep='\t', header=TRUE)
load(Rdata)

########################################
### Weighted Fraction of Reads (WFR) ###
########################################

# remove cytosines that are 0 for all
ind<-rowSums(data.cov)
dat<-data.cov[-which(ind == 0),]
idat<-tngbs.info[-which(ind == 0),]

#summarize for each assay
sum.dat<-aggregate(dat, by=list(Category=idat$V1), FUN=sum)
a1<-sum.dat$Category

# get assay names
assay<-as.data.frame(as.matrix(assay))
row.names(assay)<-assay$FASTA_name
assay<-assay[as.character(a1),]

# order by assay size
row.names(sum.dat)<-a1
sum.dat<-sum.dat[row.names(assay),]
sum.dat<-sum.dat[,-1]

# weight each assay by total sequencing reads observed in population
sum.dat<-data.matrix(sum.dat)
weights<-sum(sum.dat)/rowSums(sum.dat)
sum.dat.weighted<-apply(sum.dat,2,function(x) x*weights)

# plot boxplot of all assays coverage
pdat<-sum.dat[order(rowSums(sum.dat)),]
a.name<-assay[row.names(pdat),]
a.name<-a.name$Gene

pdat[pdat==0]<-1

pdf(file="Coverage_tngbs_populationbyassay.pdf")
seq1<-c(seq(1,10,10/1) %o% 10^(1:10))
seq1<-seq1[seq1 < max(pdat)]
boxplot(t(pdat), use.cols= TRUE, ylab="" , type = "p", col="grey", boxwex=0.4 , main="",
        log="y",axes=FALSE, xlab= "",pch=16, outcex=0.5,col.axis = "grey", cex = 1,
        font = 2, cex=0.5)
axis(2, at = seq1, labels = seq1)
axis(1, at = c(1:nrow(pdat)), labels = a.name,cex.axis=0.5,las=2)
mtext("Assay", side = 1, line = 3.5, cex = 1, font = 1)
mtext("Coverage", side = 2, line = 2.5, cex = 1, font = 1)
dev.off()

# weight each assay by total sequencing reads observed in population 1
sum.dat.weighted<-apply(sum.dat,2,function(x) x*as.numeric(as.character(assay$weights)))

# sum of weighted fraction >= 253bp: < 170bp
ab.200<-sum.dat.weighted[which(as.numeric(as.character(assay$width)) == 253),]
be.200<-colSums(sum.dat.weighted[which(as.numeric(as.character(assay$width)) < 170),])
cf.c.ratio<-ab.200/be.200

# barplot of WFR
dat_p<-cf.c.ratio
dat_p<-dat_p[order(dat_p)]
pdf(file="wfr_estimates.pdf")
op <- par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), 
          cex.lab = 1.5 , font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)
plot(dat_p, col = "black", pch = 21, bg = "grey", cex = 1, 
     ylim = c(0,max(dat_p)), ylab = "", xlab = "", axes = FALSE)
# median ccfDNA tubes
abline(h=0.0007521304, lty=2, col="blue")
mtext("ccfDNA", side = 4, line = 1, cex = 0.8,at=0.0007521304,las=2,col="blue")
#median early-processed edta
abline(h=0.01004979, lty=2, col="grey")
mtext("EDTA", side = 4, line = 1, cex = 0.8,at=0.01004979,las=2,col="grey")
#median serums
abline(h=0.0157225, lty=2, col="grey30")
mtext("Serum", side = 4, line = 1, cex = 0.8,at=0.015722,las=2,col="grey30")
#median PBMC
abline(h=0.0587956, lty=2, col="red")
mtext("PBMC", side = 4, line = 1, cex = 0.8,at=0.0587956,las=2,col="red")

abline(h= median(dat_p), col="black")
axis(1)
axis(2) 
par(las = 0)
mtext("Samples", side = 1, line = 2.5, cex = 1.5)
mtext("WFR", side = 2, line = 4.7, cex = 1.5)
dev.off()

write.table(cf.c.ratio,file="WRF_estimates.txt",sep='\t')






