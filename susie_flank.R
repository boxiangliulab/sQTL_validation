#! /usr/bin/Rscript
#====================================================================
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(VariantAnnotation))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(susieR))
suppressPackageStartupMessages(library(Rfast))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(GenomeInfoDb))
suppressPackageStartupMessages(library(parallel))
#====================================================================

option_list <- list(
  make_option(c("-v", "--vcf"), type="character", default= NULL,
              help="genotype vcf file (hg38)"),
  make_option(c("-b", "--bed"), type="character", default= NULL,
              help="bed files of QTLtools phenotype matrix [chr, start, end, ID, ID, strand]"),
  make_option(c("-n", "--name"), type="character", default= "test",
              help="name of output files [default %default]"),
  make_option("--flank", type="numeric", default= 500000,
              help="flank region around leader snps [default %default]"),
  make_option("--coverage", type="numeric", default= 0.95,
              help="susie coverage (confidence) [default %default]"),
  make_option(c("-c", "--covariance"), type="character", default= NULL,
              help="the covariance file of phenotypes [row.names=factor, col.names=sampleID]"),
  make_option("--normal", default= FALSE,  action = "store_true",
              help="A logical flag parameter [default %default]"),
  make_option("--rss", default= FALSE,  action = "store_true",
              help="choose susie_rss function [default %default]"),
  make_option("--L2", default= FALSE,  action = "store_true",
              help="choose L2 norm function [default %default]"),
  make_option(c("-q", "--qtl"), type="character", default= NULL,
              help="leader snps results file from QTLtools conditional pass"),
  make_option("--threads", type="numeric", default= 1,
              help="number of CPU threads [default %default]"),
  make_option(c("-o", "--output"), type="character", default= "./", 
              help="To specify output path [default %default]"))

parser <- OptionParser(usage = "%prog [options]", option_list=option_list,
                       description = "\nThe script of fine mapping with target region by susieR.",
                       epilogue="Yuhao Dong (yuhaodong@u.nus.edu), 2024.\n"
)

arguments <- parse_args(parser)

# 检查必需的参数是否提供
if (is.null(arguments$vcf) || is.null(arguments$bed)) {
  print_help(parser)
  stop("The arguments -v/--vcf and -b/--bed are required.", call. = FALSE)
}
# 检查参数依赖性
if (arguments$rss && is.null(arguments$qtl)) {
  stop("Error: --qtl must be specified when --rss is used.")
}

vcffile=arguments$vcf
bedfile=arguments$bed
flank=arguments$flank
confidence=arguments$coverage
covfile=arguments$covariance
qtlfile=arguments$qtl
path=arguments$output
output=paste0(arguments$output, arguments$name)

if (!file.exists(path)) {
  dir.create(path, recursive = TRUE)
  cat("path creat:", path, "\n")
}

message("vcffile: ", vcffile)
message("bedfile: ", bedfile)
message("flank: ", flank)
message("confidence: ", confidence)
if(length(covfile) == 1){message("covariance: ", covfile)}
message("susie_rss: ", arguments$rss)
if(arguments$rss){message("qtlfile: ", qtlfile)}
message("normal: ", arguments$normal)
message("output: ", output)

#
options(stringsAsFactors = FALSE)

print(paste0(arguments$name, " start"))
### 定义函数====================================================================
# get numeric genotype
seqinfo <- Seqinfo(genome = "hg38")
get_geno <- function(chr, pos, flank = 500000) {
  pos <- as.numeric(pos)
  snp <- GRanges(chr, IRanges(pos, pos), seqinfo = seqinfo)
  region <- trim(flank(snp, flank, both = TRUE))
  
  # 使用 ScanVcfParam 只读取需要的列
  vcf <- readVcf(vcffile, "hg38", param = ScanVcfParam(which = region))
  
  # 提取SNP信息和基因型数据
  snp_info <- rowRanges(vcf)
  geno <- setDT(as.data.frame(geno(vcf)$GT))
  
  # 保存行名（SNP ID）
  #snp_ids <- row.names(geno)
  snp_ids <- names(rowRanges(vcf))
  
  # 使用 data.table 的方式处理基因型
  geno_num <- geno[, lapply(.SD, function(x) {
    x <- as.character(x)
    missing_vals <- c(".", ".|.", "./.")
    valid_genos <- x[!x %in% missing_vals]
    most_common <- names(which.max(table(valid_genos)))
    x[x %in% missing_vals] <- most_common
    x <- gsub("\\|", "/", x)
    
    sapply(x, function(gt) {
      alleles <- as.numeric(strsplit(gt, "/")[[1]])
      if (all(alleles %in% c(0,1))) {
        return(sum(alleles))
      } else {
        return(sum(alleles != 0))
      }
    })
  })]
  
  # 转换为矩阵并设置行名
  result <- as.matrix(geno_num)
  row.names(result) <- snp_ids

  
  # 清理内存
  rm(geno, geno_num); gc()
  
  return(t(result))
}

# inverse normal transform
qtltools_normal_transform <- function(x) {
  # 移除NA值
  x_clean <- x[!is.na(x)]
  
  # 计算秩，使用"first"方法处理tied values
  r <- rank(x_clean, ties.method = "first")
  
  # 转换秩到[0,1]区间
  r_scaled <- r / (length(x_clean) + 1)
  
  # 应用逆正态变换
  x_transformed <- qnorm(r_scaled)
  
  # 将转换后的值放回原始向量
  x[!is.na(x)] <- x_transformed
  
  return(x)
}
### 读取文件====================================================================
# 读取表型数据文件
pheno <- fread(bedfile)
setnames(pheno, 1, "chr")
pheno <- pheno[chr != "chrY"]

# 读取cov文件
if(length(covfile) == 1){
  cov <- read.table(covfile, header=T, row.names=1, sep="\t")
  colnames(cov) <- gsub("[.]", "-", colnames(cov))
  cov[1,] <- as.numeric(factor(cov[1,], levels = c("male", "female")))
}

# 读取leader snps文件
qtl <- fread(qtlfile)

ids <- unique(paste0(qtl$V8, "|",qtl$V1))


results_list <- list()

### 封装主函数=================================================================
process_phenotype <- function(id, vcffile, bedfile, flank, confidence, covfile, pheno, cov, arguments) {
    phenotype <- str_split_fixed(id, "[|]", n=2)[,2]
    snp <- str_split_fixed(id, "[|]", n=2)[,1]
    
    print(paste0("Processing ", id))
    
  # 提取信息
  chr <- str_split_fixed(snp, "[:]", n=3)[,1]
  pos <- str_split_fixed(snp, "[:]", n=3)[,2]

  # 提取当前flank的基因型数据
  geno_numeric <- get_geno(chr, pos, flank)
  if(ncol(geno_numeric) == 0){return(NULL)}
  if(ncol(geno_numeric) == 1) {
        return(data.frame(
            variant_id = colnames(geno_numeric),
            pip = NA,
            cs_id = "one",
            min_abs_corr = NA,
            mean_abs_corr = NA,
            median_abs_corr = NA,
            id = phenotype,
            stringsAsFactors = FALSE
        ))
    }
  
  phe_geno_mat <- as.matrix(geno_numeric)
  mode(phe_geno_mat) <- "double"

  # 提取当前phenotype的表型数据
  phe_qtl_res <- subset(pheno, ID == phenotype)
  phe_y <- phe_qtl_res[, -c(1:6)] %>% t() %>% as.data.frame()
  row.names(phe_y) <- gsub("[.]", "-", row.names(phe_y))
  phe_y_vec <- as.numeric(phe_y[,1])
  names(phe_y_vec) <- row.names(phe_y)

  # regress out covariance
  if(length(covfile) == 1){
    phe_y_vec <- phe_y_vec[names(phe_y_vec) %in% colnames(cov)]
    cov <- cov[,colnames(cov) %in% names(phe_y_vec)]
    cov_matrix <- apply(cov, c(1,2), as.numeric)

    model <- lm(phe_y_vec ~ t(cov_matrix))
      
    phe_y_vec <- residuals(model)

  }

  # 对应样本数量
  phe_y_vec <- phe_y_vec[names(phe_y_vec) %in% row.names(phe_geno_mat)]
    
  phe_geno_mat <- phe_geno_mat[row.names(phe_geno_mat) %in% names(phe_y_vec),]
    
  phe_geno_mat <- as.matrix(phe_geno_mat)
  mode(phe_geno_mat) <- "double"

  phe_y_vec <- phe_y_vec[row.names(phe_geno_mat)]

  # qnorm变换
  if(arguments$normal){
    phe_y_vec <- qtltools_normal_transform(phe_y_vec)
  }

  # L2 norm
  if(arguments$L2){
    phe_y_vec <- phe_y_vec / sqrt(sum(phe_y_vec^2))
  }

  # susie ===========================================================
  if(arguments$rss){
    # 构建 LD matrix
    dat <- compute_suff_stat(phe_geno_mat, phe_y_vec, standardize = FALSE)
    R <- cov2cor(dat$XtX)

    # 读取QTL文件
    qtl <- fread(qtlfile)
    colnames(qtl)[c(1,8,12,13)] <- c("phenotypeID", "snpID", "pvalue", "slope")
    qtl <- subset(qtl, phenotypeID == phenotype)
    qtl <- qtl[qtl$snpID %in% colnames(phe_geno_mat),]

    # 计算自由度(n-2-协变量数量)
    if(length(covfile) == 1){
      df <- ncol(cov)-2-nrow(cov)
      n <- ncol(cov)
    }else{
      df <- ncol(pheno)-6-2
      n <- ncol(pheno)-6
    }
    # 计算slope SE
    qtl$slope_se <- abs(qtl$slope / qt(qtl$pvalue/2, df=df, lower.tail=FALSE))

    # 运行 susie_rss 分析
    print(paste0("Running susie_rss"))
        susie_res <- susie_rss(bhat = qtl$slope, 
                              shat = qtl$slope_se, 
                              R = R, 
                              n = n, 
                              var_y = var(phe_y_vec),
                              L = 10, 
                              coverage = confidence)
  }else{
    # 运行 susie 分析
    print(paste0("Running susie"))
    res_df <- NULL  # 初始化res_df在外部作用域
    susie_res <- tryCatch({
        res <- susie(X = phe_geno_mat,
                    y = phe_y_vec,
                    L = 10,
                    coverage = confidence)
        if(is.null(res$alpha) || nrow(res$alpha) == 0) {
            # 如果没有找到可信集合，返回一个空的结果
            return(NULL)
        }
        res
    }, error = function(e) {
        message("Error in susie for phenotype ", phenotype, ": ", e$message)
        return(NULL)
    })

  }
  if(is.null(susie_res)) {
    # 如果susie分析失败，创建一个空的结果
    return(data.frame(
            variant_id = colnames(phe_geno_mat),
            pip = NA,
            cs_id = "",
            min_abs_corr = NA,
            mean_abs_corr = NA,
            median_abs_corr = NA,
            id = phenotype,
            stringsAsFactors = FALSE
        ))
 }

  # 检查susie结果=====================================================
    if(!susie_res$converged) {
        message("SuSiE did not converge for id ", phenotype)
        return(NULL)
    }
    # 提取susie结果
    pip <- susie_res$pip
    credible_sets <- susie_res$sets$cs
    
    # 将结果与SNP注释信息整合
    res_df <- data.frame(variant_id = names(pip), pip = pip, stringsAsFactors = FALSE)
    # 如果有可信的SNP集,添加它们的置信度得分
    if (length(credible_sets) > 0) {
      # 添加credible set信息
      res_df$cs_id <- vapply(seq_along(res_df$variant_id), function(i) {
        set_indices <- which(sapply(credible_sets, function(set) i %in% set))
        if (length(set_indices) > 0) {
          paste(names(credible_sets)[set_indices], collapse = ",")
        } else {
          ""
        }
      }, FUN.VALUE = character(1))
      
      # 创建一个数据框来存储所有的purity信息
      purity_df <- data.frame(
        cs_id = rep(names(credible_sets), each = 3),
        purity_type = rep(c("min_abs_corr", "mean_abs_corr", "median_abs_corr"), times = length(credible_sets)),
        purity_value = c(susie_res$sets$purity$min.abs.corr, susie_res$sets$purity$mean.abs.corr, susie_res$sets$purity$median.abs.corr),
        stringsAsFactors = FALSE
      )
      
      # 对每个SNP,添加所属credible set的purity信息
      purity_cols <- c("min_abs_corr", "mean_abs_corr", "median_abs_corr")
      for (col in purity_cols) {
        res_df[[col]] <- vapply(res_df$cs_id, function(x) {
          if (x != "") {
            paste(purity_df$purity_value[purity_df$cs_id %in% strsplit(x, ",")[[1]] & purity_df$purity_type == col], collapse = ",")
          } else {
            ""
          }
        }, FUN.VALUE = character(1))
      }
    }
    
    # 添加表型信息
    res_df$id <- phenotype
    return(res_df)
}
### 并行处理========================================================================
message("Number of threads: ", arguments$threads)

if(arguments$threads > 1) {
    # 创建集群
    cl <- makeCluster(arguments$threads)
    
    # 导出需要的包和函数到集群
    clusterEvalQ(cl, {
        library(stringr)
        library(VariantAnnotation)
        library(dplyr)
        library(data.table)
        library(susieR)
        library(Rfast)
        library(GenomicRanges)
        library(GenomeInfoDb)
    })
    
    # 导出需要的变量到集群
    clusterExport(cl, c("vcffile", "bedfile", "flank", "confidence", 
                       "covfile", "pheno", "cov", "arguments", 
                       "get_geno", "qtltools_normal_transform", 
                       "process_phenotype", "seqinfo"))
    
    # 并行执行
    results_list <- parLapply(cl, ids, function(id) {
        tryCatch({
            process_phenotype(id, vcffile, bedfile, flank, confidence, 
                            covfile, pheno, cov, arguments)
        }, error = function(e) {
            message("Error processing ", id, ": ", e$message)
            return(NULL)
        })
    })
    
    # 关闭集群
    stopCluster(cl)
    
    # 移除NULL结果
    results_list <- results_list[!sapply(results_list, is.null)]
    
} else {
    # 单线程处理
    results_list <- list()
    for(id in ids) {
        results_list[[length(results_list) + 1]] <- process_phenotype(
            id, vcffile, bedfile, flank, confidence, covfile, pheno, cov, arguments
        )
    }
}

# 合并结果并继续后续处理
results_df <- rbindlist(results_list, fill = TRUE, use.names = TRUE)
rm(results_list);gc()

tt <- subset(results_df, cs_id > 0)
non <- subset(results_df, cs_id == "")
# 写出结果
if(arguments$rss){
  fwrite(results_df, paste0(output,".susie_rss_fine_mapping.all_snps.txt"), sep="\t")
  fwrite(tt, paste0(output,".susie_rss_fine_mapping.all_set.txt"), sep="\t")
  fwrite(tt[,'variant_id',with = T], paste0(output,".susie_rss_fine_mapping.all_set.snplist.txt"), sep="\t", col.names=F)
  fwrite(non[sample(1:nrow(non),nrow(tt)),'variant_id',with = T], paste0(output,".susie_rss_fine_mapping.no_set.control_snplist.txt"), sep="\t", col.names=F)

  }else{
  fwrite(results_df, paste0(output,".susie_fine_mapping.all_snps.txt"), sep="\t")
  fwrite(tt, paste0(output,".susie_fine_mapping.all_set.txt"), sep="\t")
  fwrite(tt[,'variant_id',with = T], paste0(output,".susie_fine_mapping.all_set.snplist.txt"), sep="\t", col.names=F)
  fwrite(non[sample(1:nrow(non),nrow(tt)),'variant_id',with = T], paste0(output,".susie_fine_mapping.no_set.control_snplist.txt"), sep="\t", col.names=F)

}


print(paste0(arguments$name, " end"))





