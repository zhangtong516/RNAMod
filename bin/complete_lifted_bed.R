#!/usr/bin/env Rscript
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

lifted_bed <- args[1]
narrowPeak <- args[2]
out_file <- args[3]

d1=fread(narrowPeak,header=F)
d2=fread(lifted_bed,header=F)

d1_clean = d1[,c(4,7,8,9,10)]

out_d = merge(d2, d1_clean, by="V4") 
setcolorder(out_d, names(d1))

fwrite(out_d,file=out_file,sep="\t",quote=F,row.names=F,col.names=F)
