# scripts below aggregate the pseduoassigned reads passing SNR threshold for eaach sample

# location of mkdf.txt files
setwd("/project/RDS-FMH-DementiaCFDNA-RW/Sunny/zac_test/methylK/test/output")

# read in each tissue:interest mkdf files (assay level pseduoassigned reads to tissue:interest) & aggregate pseduoassigned reads passing threshold for eaach sample
file.names <- dir(path, pattern ="mkdf.txt")
path = getwd()
res <- list()
for (i in 1:length(file.names)){
  if (length(file.names) == 0){
  } else {
    if (length(readLines(file.names[i])) == 0){
    } else {
      df<-read.delim(file.names[i], sep = "" , header = T,
                     na.strings ="", stringsAsFactors= F,check.names = FALSE)
      # long to wide
      agg = aggregate(as.numeric(as.character(df$frac_upsdfq_thresh)),
                      by = list(df$sample_name),
                      FUN = function(x) sum(x,na.rm= T))
      colnames(agg) <- c("sample_name",gsub("_mkdf.txt","",file.names[i]))
      res[[i]] <- agg
      
    }
  }
}

df <- sapply(res, "[[", 2)
colnames(df) <- gsub("_mkdf.txt","",file.names)
row.names(df) <- res[[1]]$sample_name

# write to file
write.table(df,file="agg_mkdf.txt",sep='\t')