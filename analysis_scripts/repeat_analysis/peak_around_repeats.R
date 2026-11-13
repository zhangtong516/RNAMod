library(data.table)
library(ggpubr)
args = commandArgs(trailingOnly = TRUE)

##
repeat_file = args[1]
peak_file = args[2]
sample_id = args[3]
##
get_overlapped_peaks = function(repeat_file, peak_file, sample_id,flanking=2000) {
    repeat_dt = fread(repeat_file)
    repeat_dt[,peak_center:=as.integer((start+end)/2)]
    repeat_dt[,start_flank:=ifelse(peak_center-flanking >1, peak_center-flanking, 1)]
    repeat_dt[,end_flank:=peak_center+flanking ]

    peak_dt = fread(peak_file)
    peak_dt[,summit:=V2+V10]
    setkey(peak_dt, V1,V2,V3)

    ##

    overlapped = foverlaps(repeat_dt, peak_dt,  
        by.x = c("chr","start_flank","end_flank"),
        by.y = c("V1","V2","V3"))

    overlapped = overlapped[!is.na(V2)]
    overlapped[,repeat_id:=str_c(chr,":",start,"-",end,"_",major_class)]

    overlapped[,relativePos:=summit-peak_center] 
    overlapped[,sample_id:=sample_id]
    return(overlapped)
}

overlapped = get_overlapped_peaks(repeat_file, peak_file, sample_id)

fwrite(overlapped, str_c(sample_id, "_overlapped_peaks.tsv"))
## 
repeat_dt = fread(repeat_file)
all_repeat_count = repeat_dt[,.N, by="major_class"]
names(all_repeat_count) = c("major_class", "total_repeat_count")
setkey(all_repeat_count, major_class)

peak_count_by_class = overlapped[,length(unique(repeat_id)), by="major_class"]
names(peak_count_by_class) = c("major_class", "N")
peak_count_by_class[,sample_id:=sample_id]
peak_count_by_class[,proportion:=round(N/sum(N),3)]
peak_count_by_class[,total_repeat_count:=all_repeat_count[major_class,total_repeat_count]]
peak_count_by_class[,overlapped_repeat_proportion:=round(N/total_repeat_count,3)]
##
pdf(str_c(sample_id, "_peak_count_by_repeat_class.pdf")) 
ggbarplot(peak_count_by_class, x="major_class", y="proportion", fill="major_class") +
    ggtitle("Proportion of repeats overlapped by class")

ggbarplot(peak_count_by_class, x="major_class", y="overlapped_repeat_proportion", fill="major_class") +
    ggtitle("Proportion of repeats overlapped over all repeats")

dev.off() 
##

pdf(str_c(sample_id, "_peak_distribution_by_repeat_class.pdf")) 
# Create 50bp bins for relativePos
overlapped[, relativePos_bin := floor(relativePos / 50) * 50]
summit_distribution = overlapped[,.N, by=c("major_class","relativePos_bin")]
summit_distribution[,proportion:=N/sum(N), by="major_class"]

ggplot(summit_distribution, aes(relativePos_bin, proportion, col=major_class)) +
    geom_line(linewidth=1)  +
    xlim(-2000,2000) +
    theme_bw()
dev.off()







