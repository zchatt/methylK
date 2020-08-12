# the function reads in summarised counts of c/g from uniquely pseudoaligned reads 
# the counts from .cgcount files are summarised in long form

## ZC test
# location_cgcnt="/Users/zacc/Documents/MtSinai/5mC_BrainBlood/Paper/re_analysis/tissuecell_analysis/test"
# cg_count<-read_cgcnt(location_cgcnt)

read_cgcnt<-function(location_cgcnt){
  file.names<-list.files(location_cgcnt, pattern = "cgcount", full.names=T)
  file.names<-list.files("/project/RDS-SMS-FFbigdata-RW/Epigenetics/tngbs/zacchatterton_tNGBS060219/mk_organise/test/output", pattern = "cgcount", full.names=T)
  file.names<-file.names[-grep("_false_",file.names)]
  file.names<-file.names[-grep("_true_",file.names)]
  
  res<-c("counts","sample","assay_long","u_psdcnt","assay","sample_name","psd_cell","ctga","cg" )
  for (i in 1:length(file.names)){
    if (length(readLines(file.names[i])) == 0){
      stop(paste0("file ",i, " has no lines"), call.=FALSE)
    } else {  
      df<-read.delim(file.names[i], sep = "" , header = T , nrows = 100,
                     na.strings ="", stringsAsFactors= F)
      df<-df[,!colnames(df) %in% c("X0_C","X0_G")]
      sample<-rep(basename(file.names))[i]
      df<-cbind(df,sample)
      df<-df[,-grep("counts.",colnames(df))]
    df2<-reshape(df, direction = "long", varying = colnames(df)[grep("chr",colnames(df))], v.names = "u_psdcnt", 
              idvar = c("counts","sample"), timevar = "assay", times = colnames(df)[grep("chr",colnames(df))],
              new.row.names=c(1:(ncol(df)*10000)))
    
    df2$u_psdcnt<-as.numeric(df2$u_psdcnt)
    #df2<-df2[!is.na(df2$u_psdcnt),]
    row.names(df2)<-NULL
    
    # get assay name
    s=unname(unlist(as.data.frame(strsplit(as.character(df2$assay),"_"))[1,]))
    s2=unname(unlist(as.data.frame(strsplit(as.character(df2$assay),"_"))[2,]))
    s3=unname(unlist(as.data.frame(strsplit(as.character(df2$assay),"_"))[3,]))
    assay<-paste(s,s2,s3,sep="_")
    
    # get sample name
    sample_name=gsub("_cgcount","",df2$sample)
    # get psd_cell
    psd_cell=unname(unlist(as.data.frame(strsplit(as.character(df2$assay),"_"))[6,]))
    # get CT/GA 
    ctga=unname(unlist(as.data.frame(strsplit(as.character(df2$assay),"_"))[4,]))
    # get C/G
    cg=unname(unlist(as.data.frame(strsplit(as.character(df2$assay),"_"))[7,]))
    
    df2<-cbind(df2,assay,sample_name,psd_cell,ctga,cg)
    res<-rbind(res,df2)
    }
  }
  df3<-res[-1,]
  colnames(df3)<-c("counts","sample","assay_long","u_psdcnt","assay","sample_name","psd_cell","ctga","cg" )
  df3$u_psdcnt<-as.numeric(df3$u_psdcnt)
  return(df3)
}

read_tfcgcnt<-function(location_cgcnt){
  file.names<-list.files(location_cgcnt, pattern = "cgcount", full.names=T)
  file.names<-list.files("/project/RDS-SMS-FFbigdata-RW/Epigenetics/tngbs/zacchatterton_tNGBS060219/mk_organise/test/output", pattern = "cgcount", full.names=T)
  file.names<-file.names[grep("_false_",file.names)]
  file.names<-file.names[grep("_true_",file.names)]
  
  res<-c("counts","sample","assay_long","u_psdcnt","assay","sample_name","psd_cell","ctga","cg" )
  for (i in 1:length(file.names)){
    if (length(readLines(file.names[i])) == 0){
      stop(paste0("file ",i, " has no lines"), call.=FALSE)
    } else {  
      df<-read.delim(file.names[i], sep = "" , header = T , nrows = 100,
                     na.strings ="", stringsAsFactors= F)
      df<-df[,!colnames(df) %in% c("X0_C","X0_G")]
      sample<-rep(basename(file.names))[i]
      df<-cbind(df,sample)
      df<-df[,-grep("counts.",colnames(df))]
      df2<-reshape(df, direction = "long", varying = colnames(df)[grep("chr",colnames(df))], v.names = "u_psdcnt", 
                   idvar = c("counts","sample"), timevar = "assay", times = colnames(df)[grep("chr",colnames(df))],
                   new.row.names=c(1:(ncol(df)*10000)))
      
      df2$u_psdcnt<-as.numeric(df2$u_psdcnt)
      #df2<-df2[!is.na(df2$u_psdcnt),]
      row.names(df2)<-NULL
      
      # get assay name
      s=unname(unlist(as.data.frame(strsplit(as.character(df2$assay),"_"))[1,]))
      s2=unname(unlist(as.data.frame(strsplit(as.character(df2$assay),"_"))[2,]))
      s3=unname(unlist(as.data.frame(strsplit(as.character(df2$assay),"_"))[3,]))
      assay<-paste(s,s2,s3,sep="_")
      
      # get sample name
      sample_name=gsub("_cgcount","",df2$sample)
      # get psd_cell
      psd_cell=unname(unlist(as.data.frame(strsplit(as.character(df2$assay),"_"))[6,]))
      # get CT/GA 
      ctga=unname(unlist(as.data.frame(strsplit(as.character(df2$assay),"_"))[4,]))
      # get C/G
      cg=unname(unlist(as.data.frame(strsplit(as.character(df2$assay),"_"))[7,]))
      
      df2<-cbind(df2,assay,sample_name,psd_cell,ctga,cg)
      res<-rbind(res,df2)
    }
  }
  df3<-res[-1,]
  colnames(df3)<-c("counts","sample","assay_long","u_psdcnt","assay","sample_name","psd_cell","ctga","cg" )
  df3$u_psdcnt<-as.numeric(df3$u_psdcnt)
  return(df3)
}


