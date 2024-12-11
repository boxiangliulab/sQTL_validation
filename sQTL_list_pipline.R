## Pipeline to get sQTL list

## susie fine mapping leader 1mb
#    Rscript ~/YuhaoDong/susie_flank_1kgp_multi.R \
#        -v ~/scratch/MAGE/vcf/MAF0.05/1kgp_chr6.MAF0.05.recode.vcf.gz \
#        -b ~/scratch/MAGE/bed/new_MAGE_spliceing_chr6.bed.gz \
#        -n MAGE_QTLtools_leader_1M_chr6_part2_3 \
#        --flank 1000000 \
#        --threads 50 \
#        -c ~/scratch/MAGE/sQTL_covariates.tab.gz \
#        -q ~/scratch/MAGE/QTLtools/chr6_part2_3.txt \
#        -o ~/scratch/MAGE/susie/


# SE function
calculate_slope_se <- function(slope, pvalue, n=729) {
  df <- n
  
  # 添加一个最小值限制
  pvalue <- pmax(pvalue, .Machine$double.eps)
  
  t_value <- qt(1 - pvalue/2, df)
  se <- abs(slope / t_value)
  
  return(se)
}

# calculate LD proxy around leader snps location(Linux Environment)
system("plink --bfile chr22 --ld-snp-list leader_snps.txt --r2 --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0.8 --out ld_results")


library(data.table)
library(stringr)
library(dplyr)

for(c in 7:22){
  chr <- paste0("chr", c)


  ## load data
  # load LD matrix
  ld <- fread(paste0("../LD/",chr,"_QTLtools_leaders_LD.ld"))
  ld <- as.data.frame(ld)
  ldd <- ld[!(ld$SNP_A == ld$SNP_B),]

  # load leader snps
  leader <- read.table(paste0("../MAGE/QTLtools/condition/MAGE_",chr,".independent_top.txt"),sep="\t", header=T)
  #leader <- subset(leader, maf >= 0.05)
  leader$slope_se <- calculate_slope_se(leader$slope, leader$pval_nominal, n=729)
  leader <- leader[c(8,1,18,21,17, 12, 19)]
  colnames(leader) <- c('variant_id', 'phenotype_id', 'slope', 'slope_se', 'pval_nominal', 'rank', 'leader')
  leader$leader_snp <- NA

  # load all sig snps
  qtl <- read.table(paste0("../MAGE/QTLtools/condition/all_sig/MAGE_",chr,".conditional.txt"))
  qtl <- subset(qtl, V20 == 1)
  qtl$slope_se <- calculate_slope_se(qtl$slope, qtl$pval_nominal, n=729)
  qtl <- qtl[c(8,1,18,21,17, 12, 19)]
  colnames(qtl) <- c('variant_id', 'phenotype_id', 'slope', 'slope_se', 'pval_nominal', 'rank', 'leader')
  qtl$leader_snp <- NA



  ## Extract all snps that located in the LD regions aroung leader snps
  proxy_list <- data.frame(matrix(ncol = ncol(qtl), nrow = 0))
  colnames(proxy_list) <- colnames(qtl)
  # combine variant ids and exon ids as unique ids
  ids <- paste(leader$variant_id, leader$phenotype_id, sep="|")

  for(i in unique(ids)){
    snp <- str_split_fixed(i, "[|]", n=2)[,1]
    exon <- str_split_fixed(i, "[|]", n=2)[,2]
    
    if(snp %in% ld$SNP_A){
      print(snp)
      print(exon)
      # Get proxy snps of target snps
      proxy <- ld[ld$SNP_A == snp,]$SNP_B
      # Get proxy snps information from original qtl list
      proxy_qtl <- qtl[(qtl$variant_id %in% proxy) & (qtl$phenotype_id == exon),]
      # if the target snp is a new leader snp and not contained in the original qtl list
      if(nrow(qtl[(qtl$variant_id == snp) & (qtl$phenotype_id == exon),])==0){
        proxy_qtl <- rbind(proxy_qtl, leader[(leader$variant_id == snp) & (leader$phenotype_id == exon), colnames(leader) %in% colnames(qtl)])
      }
      # get leader snps rank information
      #proxy_qtl[(proxy_qtl$variant_id %in% proxy) & (proxy_qtl$phenotype_id == exon),]$rank <- leader[(leader$variant_id == snp) & (leader$phenotype_id == exon),]$rank
      proxy_qtl[(proxy_qtl$variant_id == snp) & (proxy_qtl$phenotype_id == exon),]$rank <- min(leader[(leader$variant_id == snp) & (leader$phenotype_id == exon),]$rank)
      
      proxy_qtl$leader_snp <- snp
      print(i)
      proxy_list <- rbind(proxy_list, proxy_qtl)
    }
  }
  # save raw proxy list
  write.table(proxy_list, paste0(chr,"_proxy_raw_list.txt"), sep="\t", row.names=F, col.names=T, quote=F)
}


###==============================================================================

for(c in 7:22){
    chr <- paste0("chr", c)
  ## Dedup sQTLs, which are located in multiple LD region
  proxy_list$id <- paste(proxy_list$variant_id, proxy_list$phenotype_id, sep="_")
  # 1. set the label column to merge leader snp's name
  merge_column <- "leader_snp"

  # 2. find the columns except leader_snp & rank
  group_columns <- setdiff(names(proxy_list), c(merge_column, "rank"))

  # 3. Using dplyr to merge duplicated rows
  result <- proxy_list %>%
    group_by_at(group_columns) %>%
    summarise(
      leader_snp = paste(unique(leader_snp), collapse = ", "),
      rank = if(all(is.na(rank))) NA_real_ else max(rank, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    ungroup()

  result <- as.data.frame(result)

  # split information from id
  result$snp_chr <- str_split_fixed(result$variant_id, "[:]", n=4)[,1]
  result$snp_pos <- as.numeric(str_split_fixed(result$variant_id, "[:]", n=4)[,2])
  result$snp_rel <- str_split_fixed(result$variant_id, "[:]", n=4)[,3]
  result$snp_alt <- str_split_fixed(result$variant_id, "[:]", n=4)[,4]

  result$intron_start <- as.numeric(str_split_fixed(result$phenotype_id, "[:]", n=5)[,2])
  result$intron_end <- as.numeric(str_split_fixed(result$phenotype_id, "[:]", n=5)[,3])
  result$cluster <- str_split_fixed(result$phenotype_id, "[:]", n=5)[,4]
  result$gene <- str_split_fixed(result$phenotype_id, "[:]", n=5)[,5]
  result$intron_distance <- result$snp_pos - result$intron_start
  write.table(result, paste0(chr, "_proxy_list.txt"), sep="\t", row.names=F, col.names=T, quote=F)


  # save proxy snps names
  vv <- as.data.frame(unique(result$variant_id))
  write.table(vv, paste0(chr,"_sQTL_proxy.snplist.txt"), sep="\t", quote=F, row.names=F, col.names=F)
}


###merge spliceai==============================================================================
# using proxy snps names to perform spliceAI
system("vcftools --gzvcf ../MAGE/vcf/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr6.recalibrated_variants_change_ID.vcf.gz  --snps ../MAGE_leader/chr6.sQTL_proxy.snplist.txt  --recode --out chr6.sQTL_proxy")
system("spliceai -I chr6.sQTL_proxy.recode.vcf -O chr6.sQTL_proxy.spliceai.vcf.gz -R ~/yuhaodong_colo/software/puffin/resources/hg38.fa -A grch38")


###==============================================================================
for(c in 1:22){
  # read the proxy QTL list
  print(c)
  x <- read.table(paste0("skip_exon/chr",c,"_proxy_list.skip_exon.txt"),sep="\t",header=T)
  x$leader <- factor(x$leader, levels=c(0, 1), labels=c("Non_leader", "Leader"))
  x$rank <- x$rank+1
  x$slope_se <- calculate_slope_se(x$slope, x$pval_nominal, n=729)
  p <- x

  # read susie fine mapping
  susie <- fread(paste0("../susie/MAGE_QTLtools_leader_1M_chr",c,".susie_fine_mapping.all_snps.txt"),sep="\t", header=T)
  susie$phenotype_id <- susie$id
  susie$id <- paste(susie$variant_id, susie$phenotype_id, sep="_")

  susie <- susie[susie$id %in% p$id,]
  susie <- susie[!duplicated(susie$id),]
  susie <- as.data.frame(susie)

  s <- susie %>% select(-variant_id, -phenotype_id)

  data <- merge(p, s, by="id", all=F)

  # read spliceai results
  ai <- read.table(paste0("../splicaai/chr",c,".sQTL_proxy.spliceai.vcf.gz"),sep="\t",header=F)[,1:9]
  colnames(ai)[c(3,4,5,8)] <- c("variant_id", "ref", "alt", "spliceai")

  # merge informations
  da <- merge(data, ai[-c(1,2,6,7,9)], by="variant_id", all=F)


  da$spliceai_gene_symbol <- str_split_fixed(da$spliceai, "[|]", n=10)[,2]
  da$Delta_score_acceptor_gain <- str_split_fixed(da$spliceai, "[|]", n=10)[,3]
  da$Delta_score_acceptor_loss <- str_split_fixed(da$spliceai, "[|]", n=10)[,4]
  da$Delta_score_donor_gain <- str_split_fixed(da$spliceai, "[|]", n=10)[,5]
  da$Delta_score_donor_loss <- str_split_fixed(da$spliceai, "[|]", n=10)[,6]
  da$Delta_position_acceptor_gain <- str_split_fixed(da$spliceai, "[|]", n=10)[,7]
  da$Delta_position_acceptor_loss <- str_split_fixed(da$spliceai, "[|]", n=10)[,8]
  da$Delta_position_donor_gain <- str_split_fixed(da$spliceai, "[|]", n=10)[,9]
  da$Delta_position_donor_loss <- str_split_fixed(da$spliceai, "[|]", n=10)[,10]


  d <- subset(da, skip_exon_distance > 0)
  d <- d[order(d$phenotype_id, d$variant_id),]

  write.table(d[-c(2,18)], paste0("chr",c,"_sQTL_proxy_list.final.txt"), sep="\t", quote=F, row.names=F, col.names=T)

}

