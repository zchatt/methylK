# scripts below calculate the bisulfite conversion efficiencies of each sample using cytosine deamination rates within "Lambda_1" & "Lambda_3" amplicons.
library(reshape2)
library(plyr)
library(dplyr)
library(plotrix)

# location of CX_report.txt files
setwd("/project/DementiaCFDNA/Sunny/zac_test/methylK/test/output/bismark_results")

################################################################################# 
## 1 # Create matrix from bismark ".CX_report.txt" and calculate BC efficiency ##
#################################################################################
path = getwd()
# read _cgcount files
file.names <- dir(path, pattern ="CX_report.txt")
res<-c("chr","pos","strand","meth","total","ch_context","c_context","file_name","sample_name" )
for (i in 1:length(file.names)){
  if (length(file.names) == 0){
  } else {
    if (length(readLines(file.names[i])) == 0){
    } else {
      df<-read.delim(file.names[i], sep = "" , header = F ,
                     na.strings ="", stringsAsFactors= F,check.names = FALSE)
      file_name<-rep(basename(file.names))[i]
      sample_name<-gsub(".CX_report.txt","",file_name)
      df<-cbind(df,file_name,sample_name)
    }
  }
res<-rbind(res,df)
}
colnames(res)<-c("chromosome","position","strand","count_methylated","count_unmethylated","C-context","trinucleotide-context","file_name","sample_name" )
res<-res[-1,]

# calculate DNA methylation %
res$count_unmethylated <-as.numeric(as.character(res$count_unmethylated ))
res$count_methylated <-as.numeric(as.character(res$count_methylated ))
res$meth_prc<-res$count_methylated/(res$count_unmethylated + res$count_methylated) *100

# calculate bisulfite conversion by aggregating lambda meth_c/unmeth_c
res2<-res[res$chr %in% c("Lambda_1","Lambda_3"),]
d1<-aggregate(count_unmethylated ~ sample_name,res2,function(x) sum(x,na.rm = TRUE))
d2<-aggregate(count_methylated ~ sample_name,res2,function(x) sum(x,na.rm = TRUE))
d1$bc_eff<-100-(d2$count_methylated/(d1$count_unmethylated + d2$count_methylated) *100)

# write to file
write.table(d1[,-2],file="bceff.txt",sep='\t')

# plot bc efficiency of cohort
d1<-read.delim(file="bceff.txt",sep='\t',header=T)

dat_p<-d1
pdf(file="bc_efficiency.pdf")
op <- par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5 , font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)
plot(dat_p$bc_eff, col = "black", pch = 21, bg = "grey", cex = 0.5, ylim = c(90,100), ylab = "", xlab = "", axes = FALSE)
#abline(h=98, col="red", lty=2)
axis(1)
axis(2) 
par(las = 0)
mtext("Samples", side = 1, line = 2.5, cex = 1.5)
mtext("Bisulfite conversion efficiency (%)", side = 2, line = 4.7, cex = 1.5)
dev.off()



