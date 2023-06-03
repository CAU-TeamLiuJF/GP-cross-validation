#!/apps/local/software/program/R-4.0.2/bin/Rscript
## 返回输出文件中某两列的相关系数

## 加载需要的程序包
# cat('Loading required packages... \n\n')
suppressPackageStartupMessages(library("getopt"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("dplyr"))

## 命令行参数
spec <- matrix(
  c("file1", "I", 1, "character", "[Required] input file",
    "file2", "i", 1, "character", "[Optional] input file",
    "col1",  "1", 1, "integer",   "[Optional] column 1/1",
    "col2",  "2", 1, "integer",   "[Optional] column 2/2",
    "by1c",  "A", 1, "character",   "[Optional] /NULL",
    "by2c",  "a", 1, "character",   "[Optional] /NULL",
    "help",  "h", 0, "logical",   "This is Help!"),
  byrow = TRUE, ncol = 5)
opt <- getopt(spec = spec)

## 检查参数
if (!is.null(opt$help) || is.null(opt$file1)) {
  cat(paste(getopt(spec = spec, usage = T), "\n"))
  quit()
}

## 读取文件
data1 <- fread(opt$file1)

if (is.null(opt$file2)) {
  ## 只有单个文件
  if (is.null(opt$col1)) opt$col1 <- 2
  if (is.null(opt$col2)) opt$col2 <- 2

  cors <- cor(data1[, opt$col1], data1[, opt$col2])
} else {
  ## 文件2
  data2 <- fread(opt$file2)

  ## 输入文件类型
  frq <- all(names(data1) %in% c("CHR", "SNP", "A1", "A2", "MAF", "NCHROBS"))
  ld <- all(names(data1) %in% c("CHR_A", "BP_A", "SNP_A", "CHR_B", "BP_B", "SNP_B", "R2"))

  if (frq) {
    ## 计算相关的列统一成"R2"
    names(data1)[5] <- names(data2)[5] <- "R2"

    ## 合并
    data12 <- inner_join(data1, data2, by = c("CHR", "SNP"), suffix = c("_1", "_2"))

    ## 指定相同的参考碱基
    trans <- data12$A1_1 == data12$A2_2 & data12$A2_1 == data12$A1_2
    matchs <- data12$A1_1 == data12$A1_2 & data12$A2_1 == data12$A2_2
    if (sum(trans, matchs) != nrow(data12)) {
      cat("There are other types of conditions, please check (line91)\n")
      quit(status = 1)
    } else if (sum(trans) > 0) {
      ## 更改碱基和基因频率
      # cat("Ensure that reference allele (", sum(trans), ") are consistent...\n")
      data12$A1_2[trans] <- data12$A1_1[trans]
      data12$A2_2[trans] <- data12$A2_1[trans]
      data12$R2_2[trans] <- 1 - data12$R2_2[trans]
    }
  } else if (ld) {
    # cat("merging two ld result files...\n")
    data12 <- inner_join(data1, data2, by = c("SNP_A", "SNP_B"), suffix = c("_1", "_2"))
  } else {
    if (is.null(opt$by1c)) {
      cat("please provide --by1c\n")
      quit(status = 1)
    }
    if (is.null(opt$col2)) opt$col2 <- opt$col1
    if (is.null(opt$by2c)) opt$by2c <- opt$by1c

    ## 命名匹配列
    bycols1 <- as.numeric(unlist(strsplit(opt$by1c, " ")))
    bycols2 <- as.numeric(unlist(strsplit(opt$by2c, " ")))
    byN <- paste0("by", cols)
    names(data1)[opt$bycols1] <- names(data2)[opt$bycols2] <- byN

    ## 命名数据列
    names(data1)[opt$col1] <- names(data2)[opt$col2] <- "R2"

    data12 <- inner_join(data1, data2, by = byN, suffix = c("_1", "_2"))
  }

  ## 计算相关
  cors <- cor(data12$R2_1, data12$R2_2)
}

## 报告整体相关去情况
cat(cors, "\n")

## debug
opt <- list()
opt$file1 <- "/public/home/liujf/liwn/mbGS/QMSim/S0620/rand/breedBsq.ld"
opt$file2 <- "/public/home/liujf/liwn/mbGS/QMSim/S0620/rand/breedAsq.ld"
# opt$col1 = 5
# opt$by1c = 2