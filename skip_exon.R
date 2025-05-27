library(stringr)
library(data.table)
library(dplyr)

# find associated exon sites for sQTLs====================================================
for(i in 1:22){
    print(paste0("chr",i))
    #t <- fread("../chr22_proxy_list.txt")
    t <- read.table(paste0("/home/svu/e1327989/yuhaodong_colo/1kgp/MAGE_QTLtools/QTLtools/condition/all_sig/MAGE_chr",i,".conditional.txt"))
    t <- subset(t, V20 ==1)
    #t <- fread(paste0("../chr",i,"_proxy_list.txt"))
    colnames(t)[c(1,2,8,18)] <- c("intronID", "variantChrom", "variant_kgpID", "slope")

    t$intron_start <- str_split_fixed(t$intronID, "[:]", n=5)[,2]
    t$intron_end <- str_split_fixed(t$intronID, "[:]", n=5)[,3]
    t$cluster_id <- str_split_fixed(t$intronID, "[:]", n=5)[,4]
    t$intron_id <- paste(t$variantChrom, t$intron_start, t$intron_end, sep="_")
    t$intron_length <- abs(as.numeric(t$intron_end) - as.numeric(t$intron_start))

    ids <- data.frame(snp=character(), intronID=character(), site=character())
    for(snp in unique(t$variant_kgpID)){
        ts <- subset(t, variant_kgpID == snp)
        if(nrow(ts) < 2){next}
        if(all(ts$slope < 0) | all(ts$slope > 0)){next}
        
        for(cluster in unique(ts$cluster_id)){
            tt <- subset(ts, cluster_id == cluster)
            tt <- tt[!duplicated(tt$intron_id),]
            if(nrow(tt) < 2){next}
            if(all(tt$slope < 0) | all(tt$slope > 0)){next}

            for(start in unique(tt$intron_start)){
                sub_start <- tt[tt$intron_start == start, ]
                if(nrow(sub_start) < 2){next}
                if(all(sub_start$slope < 0) | all(sub_start$slope > 0)){next}
                
                start_slope <- sub_start[which.max(sub_start$intron_length),]$slope
                if(start_slope < 0 ){
                    st <- subset(sub_start, slope > 0)
                }else{
                    st <- subset(sub_start, slope < 0)
                }
                id <- paste(st$variantChrom, st$intron_end, sep="_")
                id <- unique(id)
                ids <- rbind(ids, data.frame(snp=snp, intronID=st$intronID, site=id))
            }

            for(end in unique(tt$intron_end)){
                sub_end <- tt[tt$intron_end == end, ]
                if(nrow(sub_end) < 2){next}
                if(all(sub_end$slope < 0) | all(sub_end$slope > 0)){next}
                
                end_slope <- sub_end[which.max(sub_end$intron_length),]$slope
                if(end_slope < 0 ){
                    et <- subset(sub_end, slope > 0)
                }else{
                    et <- subset(sub_end, slope < 0)
                }
                id <- paste(et$variantChrom, et$intron_start, sep="_")
                id <- unique(id)
                ids <- rbind(ids, data.frame(snp=snp, intronID=et$intronID, site=id))
            }


        }

    }

    ids$chr <- str_split_fixed(ids$site, "[_]", n=2)[,1]
    ids$start <- str_split_fixed(ids$site, "[_]", n=2)[,2]

    skip_bed <- ids[c(4,5,5,1,2,3)]
    skip_bed$id <- paste(skip_bed$snp, skip_bed$intronID, skip_bed$site, sep="|")
    skip_bed <- skip_bed[!duplicated(skip_bed$id),]

    write.table(skip_bed, paste0("skip_exon_sQTL_chr",i,".bed"), sep="\t", quote=F, row.names=F, col.names=F)

}


# overlap exons========================================================================
for(j in 1:22){
    system(paste0("bedtools intersect -a skip_exon_sQTL_chr",j,".bed  -b /data/shared/reference/hg38/gencode.v38.exon.bed  -wa -wb > skip_exon_sQTL_chr",j,"_results.bed"))
}





# 剪接调控位点与外显子匹配优化脚本
# 调控位点是单个位置（start=end），寻找最匹配的外显子
# 优先级：1) 外显子start或end等于调控位点 2) 端点距离最近 3) 中点距离最近

# 添加必要的库
library(data.table)
library(dplyr)
library(stringr)  # 用于字符串处理

# 函数：处理单个染色体的结果
process_chromosome_results <- function(chr_num) {
  
  # 读取bedtools intersect的结果
  result_file <- paste0("skip_exon_sQTL_chr", chr_num, "_results.bed")
  
  if (!file.exists(result_file)) {
    cat("Warning: File", result_file, "does not exist\n")
    return(NULL)
  }
  
  # 读取数据，假设bedtools intersect -wa -wb的输出格式
  df <- fread(result_file, header = FALSE, sep="\t")
  
  # 动态确定列数
  n_cols <- ncol(df)
  site_cols <- (n_cols / 2)  # 假设前一半是调控位点，后一半是外显子
  
  # 重命名列 - 确保包含足够的列名
  site_col_names <- paste0("site_", c("chr", "start", "end", 
                                      if(site_cols > 3) paste0("col", 4:site_cols) else NULL))
  exon_col_names <- paste0("exon_", c("chr", "start", "end", 
                                      if(site_cols > 3) paste0("col", 4:site_cols) else NULL))
  
  colnames(df) <- c(site_col_names, exon_col_names)
  
  # 创建唯一的调控位点ID
  df$site_id <- paste(df$site_chr, df$site_start, df$site_end, sep = "_")
  
  # 由于调控位点start=end，取site_start作为调控位点位置
  df$site_pos <- df$site_start
  
  # 提取基因ID信息
  df <- df %>%
    mutate(
      # 从site_col5中提取基因ID
      site_gene_id = str_split_fixed(str_split_fixed(site_col5, ":", n=5)[,5], "\\.", n=2)[,1],
      
      # 从exon_col5中提取基因ID（假设外显子也有类似的基因信息列）
      exon_gene_id = if("exon_col5" %in% colnames(.)) {
        str_split_fixed(str_split_fixed(str_split_fixed(exon_col7, "[;]", n=9)[,1], "[ ]", n=2)[,2], "[.]", n=2)[,1]
      } else {
        # 如果外显子没有col5，尝试从其他列获取基因信息
        # 这里需要根据实际的gencode exon bed格式调整
        NA_character_
      },
      
      # 判断是否同一个基因
      same_gene = !is.na(site_gene_id) & !is.na(exon_gene_id) & (site_gene_id == exon_gene_id)
    ) %>%
    mutate(
      # 判断是否端点完全匹配（外显子start或end等于调控位点位置）
      exact_boundary_match = (exon_start == site_pos | exon_end == site_pos),
      
      # 计算到外显子端点的最小距离
      dist_to_start = abs(exon_start - site_pos),
      dist_to_end = abs(exon_end - site_pos),
      min_boundary_distance = pmin(dist_to_start, dist_to_end),
      
      # 计算到外显子中点的距离（作为第四优先级）
      exon_center = (exon_start + exon_end) / 2,
      center_distance = abs(site_pos - exon_center),
      
      # 标记最近端点类型
      closest_boundary = ifelse(dist_to_start <= dist_to_end, "start", "end"),
      closest_boundary_pos = ifelse(closest_boundary == "start", exon_start, exon_end)
    )
  
  # 定义距离档位（可以根据需要调整阈值）
  df <- df %>%
    mutate(
      distance_tier = case_when(
        exact_boundary_match ~ 1,           # 第一档：完全匹配
        min_boundary_distance <= 100 ~ 2,  # 第二档：近距离（可调整阈值）
        TRUE ~ 3                           # 第三档：远距离
      )
    )
  
  # 按调控位点分组，应用三档优先级规则
  result <- df %>%
    group_by(site_col7) %>%
    arrange(distance_tier,                 # 首先按档位排序（1=完全匹配，2=近距离，3=远距离）
            desc(same_gene),               # 同档位内优先同基因
            min_boundary_distance,         # 同档位同基因内按端点距离排序
            center_distance) %>%           # 最后按中点距离排序
    slice(1) %>%                          # 每个调控位点只保留最优匹配
    ungroup()
  
  return(result)
}

# 批量处理所有染色体
all_results <- list()
for (j in 1:22) {
  cat("Processing chromosome", j, "...\n")
  chr_result <- process_chromosome_results(j)
  if (!is.null(chr_result)) {
    all_results[[j]] <- chr_result
  }
}

# 合并所有结果
final_result <- rbindlist(all_results, fill = TRUE)

# 添加匹配质量分类
final_result <- final_result %>%
  mutate(
    match_type = case_when(
      distance_tier == 1 & same_gene ~ "Tier1_Same_Gene_Exact",
      distance_tier == 1 & !same_gene ~ "Tier1_Different_Gene_Exact", 
      distance_tier == 2 & same_gene ~ "Tier2_Same_Gene_Close",
      distance_tier == 2 & !same_gene ~ "Tier2_Different_Gene_Close",
      distance_tier == 3 & same_gene ~ "Tier3_Same_Gene_Distant",
      distance_tier == 3 & !same_gene ~ "Tier3_Different_Gene_Distant",
      TRUE ~ "Unknown"
    )
  )

final_result$exon_id <- str_split_fixed(str_split_fixed(final_result$exon_col7, "[;]", n=9)[,8], "[ ]", n=3)[,3]

# 输出结果
write.table(final_result, "optimized_splicing_sites_exon_matches.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)


match_stats <- table(final_result$match_type)
print(match_stats)


# extract exon information from results================================================


final_result <- fread("skip_exon/optimized_splicing_sites_exon_matches.txt", sep="\t", header=T)

for(j in 1:22){
    #x <- fread(paste0("skip_exon/skip_exon_sQTL_chr",j,"_results.bed"),sep="\t")
    #x <- as.data.frame(x)
    #x$gene_id1 <- str_split_fixed(str_split_fixed(x$V5, "[:]", n=5)[,5], "[.]", n=2)[,1]
    #x$gene_id2 <- str_split_fixed(str_split_fixed(str_split_fixed(x$V14, "[;]", n=9)[,1], "[ ]", n=2)[,2], "[.]", n=2)[,1]
    #x$exon_id <- str_split_fixed(str_split_fixed(x$V14, "[;]", n=9)[,8], "[ ]", n=3)[,3]
    x <- subset(final_result, site_chr == paste0("chr",j))
    xx <- x[x$site_gene_id == x$exon_gene_id | x$distance_tier == 1,]
    xx <- xx[!duplicated(xx$site_col7),]


    xx$snp_pos <- str_split_fixed(xx$site_col4, "[:]", n=4)[,2]
    xx$exon_distance <- abs(as.numeric(xx$snp_pos) - as.numeric(xx$site_start))

    colnames(xx)[c(4,5,6,9,10,12,30,31)] <- c("variant_id","phenotype_id","associated_exon_site","skip_exon_start","skip_exon_end","strand","skip_exon_id","skip_exon_distance")
    
    xx <- as.data.frame(xx)[c(4,5,6,9,10,12,30,31,29)]

    xx$id <- paste(xx$variant_id, xx$phenotype_id, sep="_")

    t <- fread(paste0("../proxy_list/chr",j,"_proxy_list.txt"))

    t <- as.data.frame(t)


    data <- merge(t, xx[-c(1:2)], by="id", all=F)

    write.table(data, paste0("update_chr",j,"_proxy_list.skip_exon.txt"),sep="\t", quote=F, col.names=T, row.names=F)

}



###
da <- data[1,]
for(j in 1:22){
    data <- read.table(paste0("update_chr",j,"_proxy_list.skip_exon.txt"), header=T, sep="\t")
    da <- rbind(da, data)
}
da <- da[-1,]

write.table(da, "update_all_proxy_list.skip_exon.txt",sep="\t", quote=F, col.names=T, row.names=F)


















for(snp in unique(t$variant_kgpID)){
    ts <- subset(t, variant_kgpID == snp)
    if(length(unique(ts$intron_id)) > 10 & !(all(ts$slope > 0)) & !(all(ts$slope < 0))){
        break
    }
}

