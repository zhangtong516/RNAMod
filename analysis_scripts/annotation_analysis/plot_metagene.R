# 加载所需包
library(rtracklayer)
library(stringr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)  # 替换为您的基因组版本，如 mm10
library(ggplot2)
library(MetaTX)

args=commandArgs(trailingOnly = TRUE) 

mod_type = args[1] 
peak_dir = args[2]  ## dir to the lifted over bed files.
  
peak_files = Sys.glob(str_c(peak_dir, "/*", mod_type, "*lifted.bed"))
sample_names = unlist(lapply(str_split(basename(peak_files),paste0("__",mod_type)), function(x) x[1]))
peak_list = as.list(peak_files)
names(peak_list) = sample_names 

print(peak_list)

##
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

get_summit_gr = function(narrowpeak_file) {
    peaks <- import(narrowpeak_file, format = "narrowPeak")  # 使用 narrowPeak 格式
    summits <- GRanges(
        seqnames = seqnames(peaks),
        ranges = IRanges(
            start = start(peaks) + mcols(peaks)$peak,  # 峰峰位置 = start + peak offset (0-based)
            end = start(peaks) + mcols(peaks)$peak     # 1 bp 宽
        ),
        strand = strand(peaks),
        score = mcols(peaks)$score  # 可选：使用 score 作为权重
    )
    return(summits)
}
##
run_metatx_single=function(peak_list) {
    # get summit gr
    summit_list = lapply(peak_list, get_summit_gr)
    names(summit_list) = names(peak_list)
    # map sites
    remapresults_list <- remapCoordMulti(summit_list, txdb = txdb)
    names(remapresults_list) = names(summit_list)
    # plot
    all_density = rbindlist(lapply(names(remapresults_list), function(sample){
        remapped_result  = remapresults_list[[sample]]
        cur_density = metaTXplot(remapped_result,returnPlot = F)
        cur_density = data.table(cur_density) 
        cur_density[, sample := sample] 
        return(cur_density)
    }))
    return(all_density)
}

density_results = run_metatx_single(peak_list) 
fwrite(density_results, file = str_c(peak_dir, "/", mod_type, "_density_results.csv")) 

## plotting parameters: 
Iter <- length(density_results[,unique(sample)]) 
geom_density_adjust = 0.15 
num_bin = 10
relativeProportion = c(1, 1, 1, 1)
lambda = 2
num_bin_fiveUTR        <- num_bin
num_bin_cds            <- num_bin
num_bin_threeUTR       <- num_bin
num_bin_promoter       <- num_bin
num_bin_tail           <- num_bin
num_bin_sum            <- num_bin * 5

coord_interval         <- relativeProportion / min(relativeProportion)

coord_promoter_interval<- coord_interval[4] * num_bin
coord_promoter_Start   <- coord_interval[4] / 2
coord_promoter_end     <- coord_promoter_interval - coord_interval[4] / 2
coord_promoter         <- seq(coord_promoter_Start, coord_promoter_end,  coord_interval[4])

coord_5utr_interval    <- coord_interval[1] * num_bin
coord_5utr_Start       <- coord_promoter_interval + coord_interval[1] / 2
coord_5utr_end         <- coord_promoter_interval + coord_5utr_interval - coord_interval[1] / 2
coord_5utr             <- seq(coord_5utr_Start, coord_5utr_end, coord_interval[1])

coord_cds_interval     <- coord_interval[2] * num_bin
coord_cds_Start        <- coord_promoter_interval + coord_5utr_interval + coord_interval[2] / 2
coord_cds_end          <- coord_promoter_interval + coord_5utr_interval + coord_cds_interval - coord_interval[2] / 2
coord_cds              <- seq(coord_cds_Start, coord_cds_end, coord_interval[2])

coord_3utr_interval    <- coord_interval[3] * num_bin
coord_3utr_Start       <- coord_promoter_interval + coord_5utr_interval + coord_cds_interval + coord_interval[3] / 2
coord_3utr_end         <- coord_promoter_interval + coord_5utr_interval + coord_cds_interval + coord_3utr_interval - coord_interval[3] / 2
coord_3utr             <- seq(coord_3utr_Start, coord_3utr_end, coord_interval[3])

coord_tail_interval    <- coord_interval[4] * num_bin
coord_tail_Start       <- coord_promoter_interval + coord_5utr_interval + coord_cds_interval + coord_3utr_interval + coord_interval[4] / 2
coord_tail_end         <- coord_promoter_interval + coord_5utr_interval + coord_cds_interval + coord_3utr_interval + coord_tail_interval - coord_interval[4] / 2
coord_tail             <- seq(coord_tail_Start, coord_tail_end, coord_interval[4])

xintercept_1           <- c(coord_promoter_interval, 
                            coord_promoter_interval + coord_5utr_interval, 
                            coord_promoter_interval + coord_5utr_interval + coord_cds_interval,
                            coord_promoter_interval + coord_5utr_interval + coord_cds_interval + coord_3utr_interval)  
coord <- c(coord_promoter, coord_5utr, coord_cds, coord_3utr, coord_tail)

width_coord            <- coord_promoter_interval + coord_5utr_interval + coord_cds_interval + coord_3utr_interval + coord_tail_interval
coord                  <- coord / width_coord
xintercept_1           <- xintercept_1 / width_coord
coord_promoter_interval<- coord_promoter_interval / width_coord
coord_5utr_interval    <- coord_5utr_interval / width_coord
coord_cds_interval     <- coord_cds_interval / width_coord
coord_3utr_interval    <- coord_3utr_interval / width_coord
coord_tail_interval    <- coord_tail_interval / width_coord

for(k in 1:Iter){
    coord_interval_seq   <- c(replicate(num_bin,coord_promoter_interval / num_bin),
                                replicate(num_bin,coord_5utr_interval / num_bin),
                                replicate(num_bin,coord_cds_interval / num_bin),
                                replicate(num_bin,coord_3utr_interval / num_bin),
                                replicate(num_bin,coord_tail_interval / num_bin))
    
    max_coord            <- max(density_results$value)
    max_coord_relative   <- max_coord / sum(density_results[,value] * coord_interval_seq) 
    max_coord_relative   <- max_coord_relative / geom_density_adjust  

}


p1 <- ggplot(density_results , aes(x=coord, group=sample, weight=value)) +
    theme(panel.background =element_blank(),
        panel.grid.major = element_line(colour = 'grey', linetype = 9, size = 0.2),
        axis.text.x = element_blank(), axis.ticks = element_blank(),
        line = element_line(colour = "white", size = 0.5, linetype = 9, lineend = "butt"),
        legend.position = 'bottom', 
        axis.text=element_text(size = 8),
        legend.text=element_text(size = 8),
        axis.title.y =element_text(size = 8 ,hjust=0.5),
        title = element_text(size = 8, face='bold')) + 
    geom_density(adjust = geom_density_adjust, aes(fill=factor(sample), colour = factor(sample)),alpha = 0.15) +
    xlab("") +
    ylab("Density
        ") +
    annotate("text", x =  xintercept_1[1] / 2, y = -0.007 * max_coord_relative, label = "Promoter", size = 3) +
    annotate("text", x =  xintercept_1[1] + coord_5utr_interval / 2, y = -0.007 * max_coord_relative, label = "5'UTR", size = 3) +
    annotate("text", x =  xintercept_1[2] + coord_cds_interval / 2, y = -0.007 * max_coord_relative, label = "CDS", size = 3) +
    annotate("text", x =  xintercept_1[3] +  coord_3utr_interval / 2, y = -0.007 * max_coord_relative, label = "3'UTR", size = 3) +
    annotate("text", x =  xintercept_1[4] +  coord_tail_interval / 2, y = -0.007* max_coord_relative, label = "Tail", size = 3) + 
    geom_vline(xintercept= xintercept_1, linetype = 9, size = 0.2,colour = 'black') +
    annotate("rect", xmin = xintercept_1[1], xmax = xintercept_1[2], ymin = -0.0032 * max_coord_relative, ymax = -0.0025 * max_coord_relative, alpha = .8, colour = "black")+
    annotate("rect", xmin = xintercept_1[2], xmax = xintercept_1[3], ymin = -0.0042 * max_coord_relative, ymax = -0.0017 * max_coord_relative, alpha = .3, colour = "black")+
    annotate("rect", xmin = xintercept_1[3], xmax = xintercept_1[4], ymin = -0.0032 * max_coord_relative, ymax = -0.0025 * max_coord_relative, alpha = .8, colour = "black")

pdf(str_c(peak_dir, "/", mod_type, "_metagene_plot.pdf"), width=12, height=8) 
print(p1)
dev.off()