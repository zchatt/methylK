# R script to take the pseudocounts from kallisto and convert into dataframe for fraction of reads for both cell and tissue-type
library(data.table)

#######################
# read in data - cell #
#######################
format_psd_cell <- function(location_psd_cell) {
file.names<-list.files(location_psd_cell, pattern = ".psdccnt", full.names=T)
file.names<-file.names[-grep("total.psdccnt",file.names)]

if (all(file.exists(file.names)) == FALSE) {
  stop("Cannot find files. Please check these files have been produced", call.=FALSE)
}

res<-c("count","cell","lib")
for (i in 1:length(file.names)){
  if (length(readLines(file.names[i])) == 0){
    stop(paste0("file ",i, " has no lines"), call.=FALSE)
  } else {  
    #    N<-length(file.names)
    df<-read.delim(file.names[i], sep = "" , header = F , nrows = 100,
                   na.strings ="", stringsAsFactors= F)
    df<-cbind(df,rep(sub('.*/', '', file.names))[i])
    res<-rbind(res,df)
    # res[[i]]<-df[, -grep("counts", colnames(df))]
  }
}
df<-res[-1,]
colnames(df)<-c("count","psd_cell","lib")

# get assay name
s=unname(unlist(as.data.frame(strsplit(as.character(df$psd_cell),":"))[1,]))
assay<-gsub("_CT.*","",s)
assay<-gsub("_GA.*","",assay)

# get tissue name
s=unname(unlist(as.data.frame(strsplit(as.character(df$psd_cell),":"))[3,]))
tissue<-gsub("*methylotype_","",s)

# get sample name
s=as.character(df$lib)
sample_name<-gsub(".psdccnt","",s)

# make data.table and summarise
df<-as.data.table(cbind(sample_name,assay,tissue,df))
df$count<-as.numeric(df$count)
df<-df[, sum(count),by=list(assay,tissue,sample_name)]
colnames(df)<-c("assay","tissue","sample_name","count")
df_c<-df

#############################
# read in data - cell total #
#############################
file.names<-list.files(location_psd_cell, pattern = "_total.psdccnt", full.names=T)
if (all(file.exists(file.names)) == FALSE) {
  stop("Cannot find files. Please check these files have been produced", call.=FALSE)
}

res<-c("count","cell","lib")
for (i in 1:length(file.names)){
  if (length(readLines(file.names[i])) == 0){
    stop(paste0("file ",i, " has no lines"), call.=FALSE)
  } else {  
    #    N<-length(file.names)
    df<-read.delim(file.names[i], sep = "" , header = F , nrows = 100,
                   na.strings ="", stringsAsFactors= F)
    df<-cbind(df,rep(sub('.*/', '', file.names))[i])
    res<-rbind(res,df)
    # res[[i]]<-df[, -grep("counts", colnames(df))]
  }
}
df<-res[-1,]
colnames(df)<-c("count","psd_cell","lib")

# get assay name
s=unname(unlist(as.data.frame(strsplit(as.character(df$psd_cell),":"))[1,]))
assay<-gsub("_CT.*","",s)
assay<-gsub("_GA.*","",assay)

# get sample name
s=as.character(df$lib)
sample_name<-gsub("_total.psdccnt","",s)

# make data.table and summarise
df<-as.data.table(cbind(sample_name,assay,df))
df$count<-as.numeric(df$count)
df<-df[, sum(count),by=list(assay,sample_name)]
colnames(df)<-c("assay","sample_name","count")
df_ct<-df

# calculate fraction of reads
df_c$frac_total_pseudo<-df_c$count
p<-unique(df_c$sample_name)
for (i in 1:length(unique(p))){
  a2<-df_ct[df_ct$sample_name == p[i],]
  # for each assay 
  ua<-unique(a2$assay)
  for (z in 1:length(unique(ua))){
    c2<-df_c$count[df_c$assay == as.character(ua[z]) & df_c$sample_name == as.character(p[i])]/a2$count[a2$assay == ua[z]]
    df_c$frac_total_pseudo[df_c$assay == as.character(ua[z]) & df_c$sample_name == as.character(p[i])]<-c2
}}
# assign new name
df_cell<-df_c
return(df_cell)
}
