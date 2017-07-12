# BedAnnoGene
Given a bed file with chr and region,annotate the overlapping gene 
[usage]: ```BedAnnoGene.R bedfile gtffile outputfile```
    bedfile format: chr start end information(Arbitrary but can not be lacked)
    geffile: gtf file downloaded from GENCODE
    outputfile: file to be writen out
    needed package: data.table 1.10.4
