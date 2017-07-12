#!/bin/Rscript

library(data.table)

arg <- commandArgs(T)
if (length(arg) != 3) {
    message("[usage]: BedAnnoGene.R bedfile gtffile outputfile")
    message("    bedfile format: chr start end information(Arbitrary but can not be lacked)")
    message("    GTFfile: gtf file downloaded from GENCODE")
    message("    outputfile: file to be writen out")
    message("    needed package: data.table 1.10.4")
    stop("Please check your arguments!")
}else{
    bedfile <- arg[1]
    annofile <- arg[2]
    outfile <- arg[3]

#read file 
anno <- fread(annofile,sep="\t",header=F,skip=5)
bed <- fread(bedfile,sep="\t",header=F)
setnames(anno,c("V1","V2","V3","V4","V5","V9"),c("Chr","Gene","Type","Start","End","Info"))
anno <- anno[Type=="gene",.(Chr,Start,End,Gene=sapply(strsplit(tstrsplit(Info,";")[3][[1]],"\""),function(x)x[2]))]
setkey(anno,Chr,Start,End)
setkey(bed,V1,V2,V3)

#find overlaps by Chr
lst <- list()
for (ChrI in intersect(unique(bed$V1),unique(anno$Chr))){
  anno_reg <- anno[Chr == ChrI,.(Start,End)]
  bed_reg <- bed[V1 == ChrI,.(V2,V3)]
  setkey(anno_reg,Start,End)
  setkey(bed_reg,V2,V3)
  overl <- foverlaps(bed_reg,anno_reg,which=TRUE,nomatch = 0)
  if (nrow(overl) > 0){
    lst[[ChrI]] <- data.table(Chr=ChrI,bed[V1 == ChrI,][overl[["xid"]],.(V2,V3,V4)],anno[Chr == ChrI][overl[["yid"]],.(Gene)])
  }
}
merge_dt <- rbindlist(lst)
setnames(merge_dt,c("V2","V3","V4"),c("Start","End","Name"))

#if one region has more than one gene
torm <- list()
for (i in 1:(nrow(merge_dt)-1)){if(merge_dt[i,"Name"]==merge_dt[i+1,"Name"]){set(merge_dt,i+1L,ncol(merge_dt),paste(merge_dt[i,"Gene"],merge_dt[i+1,"Gene"],sep=";"));torm <- c(torm,list(i))}}
torm <- unlist(torm)
merge_dt <- merge_dt[-torm,]

fwrite(merge_dt,file=outfile)
}
