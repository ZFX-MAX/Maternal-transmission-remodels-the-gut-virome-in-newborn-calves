source("../R/00_config.R")
source("../R/01_utils.R")

library(data.table)
library(stringr)
meta <- fread("D:/组会/进展汇报/maternal-virome-transmission/data/meta.tsv")
library(data.table)

# 1. 筛选宝鸡数据
bj <- meta[Area == "Baoji"]

# 2. 定义 5 种配对规则 (组合名, 组1, 组2, 关联键)
cfg <- data.table(
  Comp = c("M0dF_C0dME", "M0dF_C0dOR", "M0dR_C0dME", "M0dR_C0dOR", "C0dME_C0dOR"),
  G1   = c("M0dF", "M0dF", "M0dR", "M0dR", "C0dME"),
  G2   = c("C0dME", "C0dOR", "C0dME", "C0dOR", "C0dOR"),
  Key  = c("Paired_number", "Paired_number", "Paired_number", "Paired_number", "Individuality")
)

# 3. 批量合并并提取标准列
res <- rbindlist(lapply(1:nrow(cfg), function(i) {
  d <- cfg[i]
  m <- merge(bj[Group == d$G1], bj[Group == d$G2], by = d$Key, suffixes = c("_1", "_2"))
  m[, .(Comp = d$Comp, Pair_ID = m[[d$Key]], Sample_1 = Sample_ID_1, Group_1 = Group_1, 
        Sample_2 = Sample_ID_2, Group_2 = Group_2)]
}))
fwrite(res,"D:/组会/进展汇报/maternal-virome-transmission/data/baoji_pairs.txt",sep = "\t",quote = F)
