library(data.table)
library(stringr)
library(ggpubr)

##
args = commandArgs(trailingOnly = TRUE)
peak_dir = args[1]
mod=args[2] ## m6A or m7G
analysis_id = paste0(mod,"_peaks")

overlapped_peaks_files = Sys.glob(str_c(peak_dir, "/*", mod, "*_overlapped_peaks.tsv"))

all_overlapped = rbindlist(lapply(overlapped_peaks_files, function(cur_file) {
    d=fread(cur_file)
}))

###
peak_count_by_class = all_overlapped[,length(unique(repeat_id)), by=c("sample_id", "major_class")]
names(peak_count_by_class) = c("sample_id", "major_class", "N")
peak_count_by_class[,proportion:=N/sum(N), by="sample_id"] 
# peak_count_by_class_wide = dcast(peak_count_by_class, sample_id~major_class, value.var="proportion")
# peak_count_by_class_for_radar = as.data.frame(peak_count_by_class_wide[,-"sample_id"])
# rownames(peak_count_by_class_for_radar) = peak_count_by_class_wide$sample_id
peak_count_by_class[,labels:=ifelse(proportion>=0.08, N, "")]

pdf(str_c(peak_dir, "/", analysis_id, "_peak_count_by_class.pdf"), width=12,height=8)
p1=ggbarplot(peak_count_by_class, x="sample_id", y="N",fill="major_class", 
    label = peak_count_by_class$labels, lab.pos="in", lab.size=1.5) +
    theme_bw() +
    rotate_x_text(90) 
p2=ggbarplot(peak_count_by_class, x="sample_id", y="proportion",fill="major_class")+
    theme_bw() +
    rotate_x_text(90) 

ggarrange(p1,p2, ncol=1, common.legend = TRUE, align = "v")
ggbarplot(peak_count_by_class, x="major_class", y="proportion",
        fill="sample_id",position = position_dodge()) +
    theme_bw() +
    rotate_x_text(90) 
dev.off() 

## similarity:
peak_count_by_class_wide = dcast(peak_count_by_class, sample_id~major_class, value.var="proportion")
peak_count_by_class_wide_df = peak_count_by_class_wide[,-"sample_id"]
rownames(peak_count_by_class_wide_df) = peak_count_by_class_wide$sample_id 


pheatmap::pheatmap(peak_count_by_class_wide_df, 
    file = str_c(peak_dir, "/", analysis_id, "_repeats_counts_heatmap.pdf"))

peak_count_by_class_wide_df_transpose = t(peak_count_by_class_wide_df)
colnames(peak_count_by_class_wide_df_transpose) = rownames(peak_count_by_class_wide_df)
cos_m = lsa::cosine(as.matrix(peak_count_by_class_wide_df_transpose))
pheatmap::pheatmap(cos_m, 
    file = str_c(peak_dir, "/", analysis_id, "_repeats_similarity_heatmap.pdf"))

###
bin_size=100
all_overlapped[, relativePos_bin := floor(relativePos / bin_size) * bin_size]
summit_distribution = all_overlapped[,.N, by=c("major_class","relativePos_bin", "sample_id")]
summit_distribution[,proportion:=N/sum(N), by=c("major_class", "sample_id")]

# repeats_of_interest = c("LTR/ERV", "SINE", "LINE", "DNA_transposon")
# summit_distribution = summit_distribution[major_class %in% repeats_of_interest]

pdf(str_c(peak_dir, "/", analysis_id, "_peak_distribution_around_class.pdf"), width=8,height=12)

ggplot(summit_distribution, aes(relativePos_bin, proportion, col=major_class)) +
    geom_line(linewidth=.5)  +
    xlim(-2000,2000) +
    facet_wrap(sample_id ~., ncol=1) + 
    scale_fill_brewer(palette = "Set1") +
    theme_bw()
dev.off()