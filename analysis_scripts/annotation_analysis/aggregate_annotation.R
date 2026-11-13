library(data.table)
library(stringr)
library(ggpubr)

##
anno_dir = ""

anno_files = Sys.glob(str_c(anno_dir, "/*_annotated_peaks.csv"))

all_anno = rbindlist(lapply(anno_files, function(cur_anno_file) {
    sample_name = str_extract(basename(cur_anno_file), ".+(?=_annotated_peaks.csv)")
    cur_anno = fread(cur_anno_file)
    cur_anno[, sample := sample_name]
    return(cur_anno)
}))

all_anno[str_starts(annotation, "Exon"),anno_class:="Exon"]
all_anno[str_starts(annotation, "Intron"),anno_class:="Intron"]
all_anno[is.na(anno_class), anno_class:=annotation]

##
count_by_anno_class =  all_anno[,.N, by=c("sample", "anno_class")]
count_by_anno_class[,Proportion:=N/sum(N), by="sample"]