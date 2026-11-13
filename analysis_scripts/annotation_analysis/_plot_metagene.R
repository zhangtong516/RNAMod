# 加载所需包
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)  # 替换为您的基因组版本，如 mm10
library(GenomicRanges)
library(GenomicFeatures)
library(ggplot2)
library(MetaTX)

# 步骤 1: 读取 narrowPeak 文件并提取峰峰位置（summit，第 10 列）
narrow_peak_file <- "your_narrowPeak.narrowPeak"  # 替换为您的文件路径
peaks <- import(narrow_peak_file, format = "narrowPeak")  # 使用 narrowPeak 格式

# 提取峰峰位置（1 bp 宽，用于密度计算）
summits <- GRanges(
  seqnames = seqnames(peaks),
  ranges = IRanges(
    start = start(peaks) + mcols(peaks)$peak,  # 峰峰位置 = start + peak offset (0-based)
    end = start(peaks) + mcols(peaks)$peak     # 1 bp 宽
  ),
  strand = strand(peaks),
  score = mcols(peaks)$score  # 可选：使用 score 作为权重
)

# 保存峰会为 BED 文件（可选，用于验证）
summit_bed_file <- "summits.bed"
export(summits, summit_bed_file, format = "BED")

# 步骤 2: 手动准备基因组特征（5 个部分：上游、5' UTR、CDS、3' UTR、下游）
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# 提取所有转录本
all_transcripts <- transcripts(txdb, use.names=T)

# 过滤蛋白编码转录本：使用 cdsBy 识别有 CDS 的转录本
cds_list <- cdsBy(txdb, by = "tx", use.names = TRUE)  # 获取每个转录本的 CDS（GRangesList）
has_cds <- lengths(cds_list) > 0  
protein_coding_ids <- names(cds_list[has_cds])  # 提取有 CDS 的转录本 ID

transcripts <- all_transcripts[protein_coding_ids]  # 过滤转录本

# 提取 5' UTR、CDS 和 3' UTR（仅针对蛋白编码转录本）
# 使用安全的子集操作，只选择存在的转录本ID
five_utr_all <- fiveUTRsByTranscript(txdb, use.names = TRUE)
five_utr <- five_utr_all[intersect(names(five_utr_all), protein_coding_ids)]

cds_all <- cdsBy(txdb, by = "tx", use.names = TRUE)
cds <- cds_all[intersect(names(cds_all), protein_coding_ids)]

three_utr_all <- threeUTRsByTranscript(txdb, use.names = TRUE)
three_utr <- three_utr_all[intersect(names(three_utr_all), protein_coding_ids)]

# 函数：为每个转录本创建 5 个部分（上游 2000 bp、5' UTR、CDS、3' UTR、下游 2000 bp）
# 修改 create_5parts 函数以处理缺失的ID
create_5parts <- function(tx, five_utr, cds, three_utr, upstream = 2000, downstream = 2000) {
  tx_id <- names(tx)
  n_tx <- length(tx_id)
  
  # 预先计算所有转录本的坐标信息（向量化操作）
  tx_seqnames <- seqnames(tx)
  tx_starts <- start(tx)
  tx_ends <- end(tx)
  tx_strands <- strand(tx)
  
  # 向量化计算TSS和TTS位置
  is_plus <- tx_strands == "+"
  tss_pos <- ifelse(is_plus, tx_starts, tx_ends)
  tts_pos <- ifelse(is_plus, tx_ends, tx_starts)
  
  # 批量创建上游区域
  upstream_gr <- GRanges(
    seqnames = tx_seqnames,
    ranges = IRanges(start = pmax(1, tss_pos - upstream), end = tss_pos - 1),
    strand = tx_strands,
    feature = "upstream",
    tx_id = tx_id
  )
  
  # 批量创建下游区域
  downstream_gr <- GRanges(
    seqnames = tx_seqnames,
    ranges = IRanges(start = tts_pos + 1, end = tts_pos + downstream),
    strand = tx_strands,
    feature = "downstream",
    tx_id = tx_id
  )
  
  # 向量化处理UTR和CDS区域
  # 提取所有5' UTR并添加metadata
  five_utr_valid <- five_utr[intersect(names(five_utr), tx_id)]
  if (length(five_utr_valid) > 0) {
    five_utr_unlisted <- unlist(five_utr_valid)
    mcols(five_utr_unlisted)$feature <- "5utr"
    mcols(five_utr_unlisted)$tx_id <- rep(names(five_utr_valid), lengths(five_utr_valid))
  } else {
    five_utr_unlisted <- GRanges()
  }
  
  # 提取所有CDS并添加metadata
  cds_valid <- cds[intersect(names(cds), tx_id)]
  if (length(cds_valid) > 0) {
    cds_unlisted <- unlist(cds_valid)
    mcols(cds_unlisted)$feature <- "cds"
    mcols(cds_unlisted)$tx_id <- rep(names(cds_valid), lengths(cds_valid))
  } else {
    cds_unlisted <- GRanges()
  }
  
  # 提取所有3' UTR并添加metadata
  three_utr_valid <- three_utr[intersect(names(three_utr), tx_id)]
  if (length(three_utr_valid) > 0) {
    three_utr_unlisted <- unlist(three_utr_valid)
    mcols(three_utr_unlisted)$feature <- "3utr"
    mcols(three_utr_unlisted)$tx_id <- rep(names(three_utr_valid), lengths(three_utr_valid))
  } else {
    three_utr_unlisted <- GRanges()
  }
  
  # 合并所有区域
  all_regions <- c(upstream_gr, five_utr_unlisted, cds_unlisted, three_utr_unlisted, downstream_gr)
  
  # 按转录本ID和特征类型排序（可选，用于保持一致的顺序）
  feature_order <- c("upstream", "5utr", "cds", "3utr", "downstream")
  all_regions$feature <- factor(all_regions$feature, levels = feature_order)
  all_regions <- all_regions[order(all_regions$tx_id, all_regions$feature)]
  
  return(all_regions)
}

# 确保所有对象都有相同的转录本ID集合
# 找到所有对象中都存在的转录本ID
common_ids <- Reduce(intersect, list(
  names(transcripts),
  names(five_utr),
  names(cds),
  names(three_utr)
))

# 使用相同的ID集合过滤所有对象
transcripts <- transcripts[common_ids]
five_utr <- five_utr[common_ids]
cds <- cds[common_ids]
three_utr <- three_utr[common_ids]

# 现在可以安全地调用 create_5parts
gf5_parts <- create_5parts(transcripts, five_utr, cds, three_utr, upstream = 2000, downstream = 2000)

# 步骤 3: 计算每个区域的归一化峰峰密度
# 函数：计算指定区域的密度（bins 数量根据区域调整）
# 优化的密度计算函数 - 向量化处理所有区域
compute_all_densities <- function(summits, regions, regions_order, bins_per_region) {
  # 预先计算所有重叠
  hits <- findOverlaps(summits, regions)
  
  # 创建结果列表
  densities <- vector("list", length(regions_order))
  names(densities) <- regions_order
  
  # 为每个区域类型批量处理
  for (i in seq_along(regions_order)) {
    region_name <- regions_order[i]
    bins <- bins_per_region[i]
    
    # 子集该类型的所有区域
    reg_indices <- which(mcols(regions)$feature == region_name)
    if (length(reg_indices) == 0) {
      densities[[i]] <- rep(0, bins)
      next
    }
    
    reg_subset <- regions[reg_indices]
    
    # 找到与该类型区域重叠的峰峰
    region_hits <- hits[subjectHits(hits) %in% reg_indices]
    
    if (length(region_hits) == 0) {
      densities[[i]] <- rep(0, bins)
      next
    }
    
    # 向量化计算相对位置和bin索引
    summit_indices <- queryHits(region_hits)
    region_indices_local <- match(subjectHits(region_hits), reg_indices)
    
    summit_starts <- start(summits)[summit_indices]
    region_starts <- start(reg_subset)[region_indices_local]
    region_widths <- width(reg_subset)[region_indices_local]
    
    # 批量计算相对位置
    rel_pos <- (summit_starts - region_starts) / region_widths
    bin_indices <- pmin(floor(rel_pos * bins) + 1, bins)
    
    # 使用tabulate快速计数
    counts <- tabulate(bin_indices, nbins = bins)
    
    # 归一化
    total_summits <- length(summits)
    avg_bin_width <- mean(width(reg_subset)) / bins
    densities[[i]] <- counts / (total_summits * avg_bin_width)
  }
  
  return(densities)
}

# 替换原来的mapply调用 - 使用50bp bins进行平滑
regions_order <- c("upstream", "5utr", "cds", "3utr", "downstream")

# 计算每个区域的平均长度，然后基于50bp计算bins数量
calculate_bins_for_50bp <- function(regions, region_name, bin_size = 50) {
  region_subset <- regions[mcols(regions)$feature == region_name]
  if (length(region_subset) == 0) return(1)
  
  avg_length <- mean(width(region_subset))
  bins <- max(1, round(avg_length / bin_size))
  return(bins)
}

# 为每个区域计算基于50bp的bins数量
bins_per_region <- sapply(regions_order, function(region_name) {
  calculate_bins_for_50bp(gf5_parts, region_name, bin_size = 50)
})

cat("Bins per region (50bp each):\n")
for (i in seq_along(regions_order)) {
  region_subset <- gf5_parts[mcols(gf5_parts)$feature == regions_order[i]]
  avg_length <- if(length(region_subset) > 0) mean(width(region_subset)) else 0
  cat(sprintf("%s: %d bins (avg length: %.1f bp)\n", 
              regions_order[i], bins_per_region[i], avg_length))
}

# 优化的密度计算函数 - 向量化处理所有区域，使用50bp bins
compute_all_densities_50bp <- function(summits, regions, regions_order, bin_size = 50) {
  # 预先计算所有重叠
  hits <- findOverlaps(summits, regions)
  
  # 创建结果列表
  densities <- vector("list", length(regions_order))
  names(densities) <- regions_order
  
  # 为每个区域类型批量处理
  for (i in seq_along(regions_order)) {
    region_name <- regions_order[i]
    
    # 子集该类型的所有区域
    reg_indices <- which(mcols(regions)$feature == region_name)
    if (length(reg_indices) == 0) {
      densities[[i]] <- numeric(0)
      next
    }
    
    reg_subset <- regions[reg_indices]
    avg_length <- mean(width(reg_subset))
    bins <- max(1, round(avg_length / bin_size))
    
    # 找到与该类型区域重叠的峰峰
    region_hits <- hits[subjectHits(hits) %in% reg_indices]
    
    if (length(region_hits) == 0) {
      densities[[i]] <- rep(0, bins)
      next
    }
    
    # 向量化计算相对位置和bin索引
    summit_indices <- queryHits(region_hits)
    region_indices_local <- match(subjectHits(region_hits), reg_indices)
    
    summit_starts <- start(summits)[summit_indices]
    region_starts <- start(reg_subset)[region_indices_local]
    region_widths <- width(reg_subset)[region_indices_local]
    
    # 批量计算相对位置
    rel_pos <- (summit_starts - region_starts) / region_widths
    bin_indices <- pmin(floor(rel_pos * bins) + 1, bins)
    
    # 使用tabulate快速计数
    counts <- tabulate(bin_indices, nbins = bins)
    
    # 归一化：每个bin代表50bp，密度为 counts per 50bp per 1000 summits
    total_summits <- length(summits)
    densities[[i]] <- counts / (total_summits / 1000) / (bin_size / 1000)  # 标准化为每1000个summits每kb的密度
  }
  
  return(densities)
}

# 使用50bp bins计算密度
densities <- compute_all_densities_50bp(summits, gf5_parts, regions_order, bin_size = 50)

# 更新bins_per_region以匹配实际计算的bins
bins_per_region <- sapply(densities, length)

# 步骤 4: 准备绘图数据（等长区域 + 分别标准化密度）
standard_length <- 1000  # 所有区域使用相同的显示长度

x_positions <- c()
density_values <- c()
region_labels <- c()

# 为每个区域分别标准化密度并设置等长显示
for (i in seq_along(regions_order)) {
  r <- regions_order[i]
  b <- bins_per_region[i]
  dens <- densities[[i]]
  
  # 分别标准化每个区域的密度（0-1标准化）
  if (max(dens) > 0) {
    normalized_dens <- (dens - min(dens)) / (max(dens) - min(dens))
  } else {
    normalized_dens <- dens
  }
  
  # 或者使用Z-score标准化
  # normalized_dens <- (dens - mean(dens)) / sd(dens)
  
  # 每个区域占用相同的显示长度
  start_pos <- (i - 1) * standard_length
  end_pos <- i * standard_length
  
  if (i == 1) {
    # 第一个区域包含起始点
    x_pos <- seq(start_pos, end_pos, length.out = b + 1)[1:b]
  } else if (i == length(regions_order)) {
    # 最后一个区域包含结束点
    x_pos <- seq(start_pos, end_pos, length.out = b + 1)
    x_pos <- x_pos[-1]  # 移除第一个点以避免重复
  } else {
    # 中间区域
    x_pos <- seq(start_pos, end_pos, length.out = b + 1)
    x_pos <- x_pos[-1]  # 移除第一个点以避免重复
  }
  
  x_positions <- c(x_positions, x_pos)
  density_values <- c(density_values, normalized_dens)
  region_labels <- c(region_labels, rep(r, length(x_pos)))
}

# 创建数据框
plot_data <- data.frame(
  Position = x_positions,
  Density = density_values,
  Region = region_labels
)

# 计算区域边界位置
boundary_positions <- (1:(length(regions_order)-1)) * standard_length

# 绘制等长区域的标准化密度图
p <- ggplot(plot_data, aes(x = Position, y = Density)) +
  geom_line(size = 1, color = "blue") +
  geom_vline(xintercept = boundary_positions, linetype = "dashed", color = "gray", alpha = 0.7) +
  
  # 为每个区域添加不同的背景色
  annotate("rect", 
           xmin = c(0, boundary_positions), 
           xmax = c(boundary_positions, length(regions_order) * standard_length),
           ymin = -Inf, ymax = Inf, 
           alpha = 0.1, 
           fill = c("red", "green", "blue", "orange", "purple")) +
  
  scale_x_continuous(
    breaks = seq(standard_length/2, length(regions_order) * standard_length - standard_length/2, by = standard_length),
    labels = c("Upstream\n(-2kb)", "5' UTR", "CDS", "3' UTR", "Downstream\n(+2kb)")
  ) +
  
  # 添加区域标签
  annotate("text", 
           x = seq(standard_length/2, length(regions_order) * standard_length - standard_length/2, by = standard_length),
           y = max(plot_data$Density) * 0.9,
           label = c("Upstream", "5' UTR", "CDS", "3' UTR", "Downstream"),
           size = 3, color = "black", fontface = "bold") +
  
  labs(title = "Normalized Metagene Plot (Equal Region Lengths)",
       subtitle = "Each region normalized separately (0-1 scale) and displayed with equal length",
       x = "Standardized Gene Features",
       y = "Normalized Peak Density (0-1 scale)") +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# 打印每个区域的原始密度统计
cat("Original density statistics by region:\n")
for (i in seq_along(regions_order)) {
  r <- regions_order[i]
  dens <- densities[[i]]
  cat(sprintf("%s: mean=%.4f, max=%.4f, min=%.4f\n", 
              r, mean(dens), max(dens), min(dens)))
}

# 可选：打印转录本数量以验证过滤
cat("Number of protein-coding transcripts:", length(transcripts), "\n")

# 替代方案：先计算密度，然后应用50bp窗口平滑
apply_smoothing_50bp <- function(density_values, original_bin_size, target_bin_size = 50) {
  if (length(density_values) == 0) return(numeric(0))
  
  # 计算平滑窗口大小
  window_size <- max(1, round(target_bin_size / original_bin_size))
  
  # 应用移动平均平滑
  if (window_size == 1) {
    return(density_values)
  }
  
  # 使用卷积进行移动平均
  kernel <- rep(1/window_size, window_size)
  
  # 处理边界：使用反射填充
  padded_values <- c(rev(density_values[1:min(window_size, length(density_values))]),
                     density_values,
                     rev(density_values[max(1, length(density_values)-window_size+1):length(density_values)]))
  
  smoothed <- stats::filter(padded_values, kernel, sides = 2)
  
  # 移除填充并返回原始长度
  start_idx <- window_size + 1
  end_idx <- start_idx + length(density_values) - 1
  
  return(as.numeric(smoothed[start_idx:end_idx]))
}

# 对每个区域的密度应用50bp平滑
smoothed_densities <- vector("list", length(regions_order))
names(smoothed_densities) <- regions_order

for (i in seq_along(regions_order)) {
  region_name <- regions_order[i]
  original_density <- densities[[i]]
  
  # 估算原始bin大小
  region_subset <- gf5_parts[mcols(gf5_parts)$feature == region_name]
  if (length(region_subset) > 0 && length(original_density) > 0) {
    avg_length <- mean(width(region_subset))
    original_bin_size <- avg_length / length(original_density)
    
    smoothed_densities[[i]] <- apply_smoothing_50bp(original_density, original_bin_size, 50)
  } else {
    smoothed_densities[[i]] <- original_density
  }
}

# 使用平滑后的密度
densities <- smoothed_densities
bins_per_region <- sapply(densities, length)