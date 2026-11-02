# ============================================================================
#  Fig1_Overview3.R
#  Fig1 全流程：病毒组 (Step 1-7) + 细菌组 (Step B1-B5)
#  统一整理版 —— 2026-09-21
# ============================================================================

# ----------------------------------------------------------------------------
# 0.1 依赖
# ----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(microeco)
  library(vegan)
  library(patchwork)
  library(RColorBrewer)
  library(matrixStats)
  library(ggpubr)          # stat_compare_means
})

# ----------------------------------------------------------------------------
# 0.2 配置与目录
# ----------------------------------------------------------------------------
source("../R/00_config.R")
source("../R/01_utils.R")

DATA_DIR <- "../data"                       # ← 统一数据目录
OUT_DIR  <- "../output/Fig1"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------------------------------------------
# 0.3 参数兜底（00_config.R 未定义时使用默认值）
# ----------------------------------------------------------------------------
default_params <- list(
  MAX_CONTAM      = 10,
  QUALITY_KEEP    = c("High-quality", "Medium-quality","Low-quality","Complete"),
  EXCLUDE_GROUP   = character(0),
  MIN_CONTIG_LEN  = 5000,
  MIN_RPKM        = 1,
  TOP_FAMILY_N    = 8,
  TOP_HOST_N      = 6,
  N_PERM          = 999,
  N_CORES         = 1,
  GROUP_LEVELS    = c("Baoji_C0dME", "Baoji_C0dOR", "Gansu_C0dME", "Gansu_C70dF",
                      "Baoji_M0dF",  "Gansu_M0dF",  "Baoji_M0dR",  "Gansu_M0dR",
                      "Baoji_M21dF", "Baoji_M21dR"),
  BAR_COLORS      = brewer.pal(8, "Set2")
)
for (nm in names(default_params)) {
  if (!exists(nm, inherits = TRUE)) assign(nm, default_params[[nm]])
}
TAX_LEVELS <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# ----------------------------------------------------------------------------
# 0.4 工具函数
# ----------------------------------------------------------------------------
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a

## 从 contig/cluster 名中提取长度: xxx_length_12345_cov_xx -> 12345
extract_contig_length <- function(x) {
  as.integer(sub(".*_length_([0-9]+)_.*", "\\1", x))
}

## 行最大值
row_max <- function(x) {
  if (requireNamespace("matrixStats", quietly = TRUE)) matrixStats::rowMaxs(x)
  else apply(x, 1, max, na.rm = TRUE)
}

## p 值格式化
fmt_p <- function(p, digits = 3) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("< 0.001")
  if (p < 0.01)  return("< 0.01")
  if (p < 0.05)  return("< 0.05")
  sprintf("= %.*f", digits, p)
}

## 出版级主题
theme_pub <- function(base_size = 12, base_family = "") {
  theme_classic(base_size = base_size, base_family = base_family) %+replace%
    theme(
      plot.title       = element_text(hjust = 0.5, face = "bold", size = rel(1.1)),
      axis.title       = element_text(size = rel(1.0)),
      axis.text        = element_text(size = rel(0.9)),
      legend.title     = element_text(face = "bold", size = rel(0.95)),
      legend.text      = element_text(size = rel(0.85)),
      strip.text       = element_text(face = "bold", size = rel(0.9)),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}

## 保存兜底
if (!exists("save_tsv")) save_tsv <- function(dt, name)
  fwrite(as.data.table(dt), file.path(OUT_DIR, name), sep = "\t")
if (!exists("save_pdf")) save_pdf <- function(p, name, w = 8, h = 6)
  ggsave(file.path(OUT_DIR, name), p, width = w, height = h, device = cairo_pdf)

## meta 与距离矩阵强制对齐
align_meta <- function(meta_dt, labels) {
  m <- as.data.table(meta_dt)
  m <- m[m$Sample_ID %in% labels]
  m <- m[match(labels, m$Sample_ID)]
  if (nrow(m) != length(labels) || any(is.na(m$Sample_ID))) {
    stop(sprintf("meta 与距离矩阵不匹配: dist=%d, matched=%d",
                 length(labels), sum(!is.na(m$Sample_ID))))
  }
  m
}

## adonis2 封装（修复版：正确求值 + 公式环境锚定）
fast_adonis2 <- function(formula, data, permutations = N_PERM,
                         parallel = 1, by = NULL, ...) {
  if (!inherits(formula, "formula")) stop("'formula' 必须是 formula 对象")
  pf <- parent.frame()
  
  ## 1. 求值左侧 → 本地变量 dm
  dm <- eval(formula[[2L]], envir = pf)
  if (!inherits(dm, "dist"))
    stop(sprintf("公式左侧必须是 dist 对象，当前是 %s", class(dm)[1]))
  
  labs <- attr(dm, "Labels")
  n    <- attr(dm, "Size")
  
  ## 2. data 对齐
  data <- as.data.table(data)
  if (!is.null(labs)) {
    if ("Sample_ID" %in% names(data)) {
      data <- align_meta(data, labs)
    } else if (!is.null(rownames(data)) && all(labs %in% rownames(data))) {
      data <- data[match(labs, rownames(data)), ]
    }
  }
  if (nrow(data) != n)
    stop(sprintf("行数不匹配: data = %d, dist = %d", nrow(data), n))
  
  ## 3. 右侧变量校验
  rhs  <- all.vars(formula[[3L]])
  miss <- setdiff(rhs, names(data))
  if (length(miss)) stop("data 中缺少右侧变量: ", paste(miss, collapse = ", "))
  for (v in rhs) {
    vv <- data[[v]]
    if (length(unique(vv[!is.na(vv)])) < 2)
      stop(sprintf("变量 %s 只有 1 个有效水平，无法做 PERMANOVA", v))
  }
  
  ## 4. 转 data.frame + 行名对齐
  dframe <- as.data.frame(data)
  rownames(dframe) <- labs %||% rownames(dframe)
  
  ## 5. 关键修复：把公式左侧替换成本地变量 dm
  ##    update(old, dm ~ .)  →  LHS 变成 dm，RHS 保持不变
  ##    然后锚定环境到本帧，vegan eval 时就能找到 dm
  formula <- update(formula, dm ~ .)
  environment(formula) <- environment()
  
  vegan::adonis2(formula, data = dframe, permutations = permutations,
                 by = by, parallel = max(1L, parallel), ...)
}
run_variance_partition <- function(dist_full, st, label = "") {
  
  dist_full <- as.matrix(dist_full)
  if (!is.null(rownames(dist_full)) && "Sample_ID" %in% names(st)) {
    st <- align_meta(as.data.table(st), rownames(dist_full))
  }
  
  ## 预先求值成具名 dist 对象（供公式引用）
  dm_all <- as.dist(dist_full)
  stopifnot(attr(dm_all, "Size") == nrow(st))
  
  ## (1) 单因素：dm 在 lapply 的匿名函数帧内创建，公式同帧 → 天然安全
  single <- rbindlist(lapply(c("Area", "Niche", "Age", "Time", "Group_Area"), function(fac) {
    if (!fac %in% names(st)) return(NULL)
    vi <- !is.na(st[[fac]])
    if (sum(vi) <= 10 || length(unique(st[[fac]][vi])) < 2) return(NULL)
    
    d  <- data.frame(x = factor(st[[fac]][vi]))
    dm <- as.dist(dist_full[vi, vi])
    r  <- vegan::adonis2(dm ~ x, data = d, permutations = N_PERM)
    
    data.table(Factor = fac, R2 = round(r$R2[1], 4),
               F_value = round(r$F[1], 3), P_value = r$`Pr(>F)`[1], N = sum(vi))
  }))
  
  ## (2) 多因素：走修复后的 fast_adonis2
  multi_vars <- intersect(c("Area", "Niche", "Age"), names(st))
  multi <- NULL
  if (length(multi_vars) >= 2) {
    fml <- as.formula(paste("dm_all ~", paste(multi_vars, collapse = " + ")))
    multi_res <- fast_adonis2(fml, data = st, by = "margin", permutations = N_PERM)
    k <- length(multi_vars)
    multi <- data.table(Factor  = rownames(multi_res)[seq_len(k)],
                        R2      = round(multi_res$R2[seq_len(k)], 4),
                        F_value = round(multi_res$F[seq_len(k)], 3),
                        P_value = multi_res$`Pr(>F)`[seq_len(k)])
  }
  
  ## (3) 分队列：同样用局部 dm
  shared <- c("C0dME", "M0dF", "M0dR")
  cohort <- rbindlist(lapply(c("Baoji", "Gansu"), function(area) {
    idx <- which(st$Group_Area %in% paste0(area, "_", shared))
    if (length(idx) <= 10) return(NULL)
    
    st_a <- copy(st[idx])
    st_a[, Group_simple := factor(Group, levels = intersect(shared, unique(Group)))]
    if (nlevels(droplevels(st_a$Group_simple)) < 2) return(NULL)
    
    dm <- as.dist(dist_full[idx, idx])
    r  <- vegan::adonis2(dm ~ Group_simple, data = st_a, permutations = N_PERM)
    
    data.table(Cohort = area, R2 = round(r$R2[1], 4),
               F_value = round(r$F[1], 3), P_value = r$`Pr(>F)`[1], N = length(idx))
  }))
  
  if (nzchar(label)) {
    message(sprintf("--- %s single-factor PERMANOVA ---", label));        print(single)
    message(sprintf("--- %s multi-factor PERMANOVA (margin) ---", label)); print(multi)
    message(sprintf("--- %s cohort-specific PERMANOVA ---", label));       print(cohort)
  }
  list(single = single, multi = multi, cohort = cohort)
}
## Alpha 多样性箱线图（带显著性标记）
plot_alpha_box <- function(alpha_dt, comparisons, colors,
                           ylab = "Shannon Diversity Index (based on RPKM)") {
  ggplot(alpha_dt, aes(Group_Area, shannon, fill = Group_Area, color = Group_Area)) +
    geom_boxplot(width = 0.6, color = "#3F3F3F", linewidth = 0.5,
                 alpha = 0.8, outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 1.2, alpha = 0.6, show.legend = FALSE) +
    facet_grid(. ~ Area, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
    labs(x = NULL, y = ylab) +
    theme_pub() +
    lapply(names(comparisons), function(area) {
      stat_compare_means(data = alpha_dt[Area == area],
                         comparisons = comparisons[[area]],
                         method = "wilcox.test", label = "p.signif",
                         inherit.aes = TRUE)
    }) +
    theme(axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
          axis.text.y  = element_text(size = 7),
          axis.title.y = element_text(size = 8),
          legend.position = "none",
          strip.text   = element_text(size = 8, face = "bold"),
          panel.grid   = element_blank())
}

## Beta 多样性一站式：dist + PERMANOVA + betadisper + PCoA
run_beta <- function(rpkm_mat, sample_table, tag = "", w = 10, h = 8,
                     make_plot = TRUE) {
  labels <- colnames(rpkm_mat)
  st     <- align_meta(sample_table, labels)
  # dm     <- vegdist(t(rpkm_mat[labels], method = "bray"))   # 占位，见下一行修正
  dm     <- vegdist(t(rpkm_mat), method = "bray")
  
  perm   <- vegan::adonis2(dm ~ Group_Area, data = st, permutations = N_PERM)
  p_permanova <- perm$`Pr(>F)`[1]
  r2_permanova <- perm$R2[1]
  p_disp <- anova(betadisper(dm, st$Group_Area))$`Pr(>F)`[1]
  
  pcoa <- cmdscale(dm, k = 2, eig = TRUE)
  ve   <- round(pcoa$eig[1:2] / sum(pcoa$eig[pcoa$eig > 0]) * 100, 2)
  pdf_ <- data.table(Sample_ID = rownames(pcoa$points),
                     PCo1 = pcoa$points[, 1], PCo2 = pcoa$points[, 2])
  pdf_ <- merge(pdf_, st[, .(Sample_ID, Group_Area, Area)], by = "Sample_ID")
  set(pdf_, NULL, "Group_Area", factor(pdf_$Group_Area, levels = GROUP_LEVELS))
  
  pl <- NULL
  if (make_plot) {
    cols <- group_color[levels(pdf_$Group_Area)]
    pl <- ggplot(pdf_, aes(PCo1, PCo2, color = Group_Area, fill = Group_Area, shape = Area)) +
      stat_ellipse(geom = "polygon", alpha = 0.15, level = 0.95, show.legend = FALSE) +
      geom_point(size = 2.5, alpha = 0.7) +
      scale_color_manual(values = cols) +
      scale_fill_manual(values = cols) +
      scale_shape_manual(values = c(16, 17)) +
      labs(x = sprintf("PCo1 (%s%%)", ve[1]),
           y = sprintf("PCo2 (%s%%)", ve[2]),
           title = sprintf("%sPCoA (PERMANOVA R\u00b2 = %.3f, P %s; Betadisper P = %.3f)",
                           tag, r2_permanova, fmt_p(p_permanova), p_disp)) +
      theme_classic() +
      theme(legend.position = "right",
            plot.title = element_text(size = 8, hjust = 0.5),
            axis.text  = element_text(size = 8),
            axis.title = element_text(size = 9))
  }
  
  sig <- data.table(Test      = c("PERMANOVA", "Betadisper"),
                    Statistic = c(sprintf("R2 = %.3f", r2_permanova),
                                  "Homogeneity of variances"),
                    P_value   = c(p_permanova, p_disp))
  
  list(dist = dm, dist_full = as.matrix(dm), st = st, perm = perm,
       pcoa_df = pdf_, var_exp = ve, plot = pl, sig = sig,
       p_permanova = p_permanova, r2_permanova = r2_permanova, p_disp = p_disp)
}

## 组间成对距离
build_pair_dist <- function(st, dist_full, g1, g2, comparison, source = NULL) {
  s1 <- intersect(st$Sample_ID[st$Group_Area == g1], rownames(dist_full))
  s2 <- intersect(st$Sample_ID[st$Group_Area == g2], rownames(dist_full))
  if (!length(s1) || !length(s2)) return(NULL)
  data.table(Comparison = comparison, Group1 = g1, Group2 = g2,
             Distance = as.vector(dist_full[s1, s2]), Source = source)
}

## 距离箱线图
DIST_COLORS <- c(M21dR = "#EDBA83", M0dR = "#FB9C1D",
                 M21dF = "#9ED5D8", M0dF = "#75CED9", Oral = "#5C9E78")

plot_dist_box <- function(dt, comparisons, title, colors = DIST_COLORS) {
  ggplot(dt, aes(Mother_Group, Distance, fill = Mother_Group, color = Mother_Group)) +
    geom_boxplot(width = 0.6, color = "gray", linewidth = 0.5,
                 alpha = 0.85, outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 1.5, alpha = 0.7, show.legend = FALSE) +
    facet_grid(. ~ Comparison, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    labs(x = "Maternal Group", y = "Bray-Curtis Distance", title = title) +
    theme_bw(base_size = 9) +
    theme(legend.position = "right",
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8),
          axis.title  = element_text(size = 9),
          strip.text  = element_text(size = 8, face = "bold"),
          plot.title  = element_text(size = 10, hjust = 0.5, face = "bold"),
          panel.grid.minor = element_blank(),
          plot.margin = margin(10, 15, 10, 10, unit = "pt")) +
    lapply(names(comparisons), function(comp) {
      stat_compare_means(data = dt[Comparison == comp],
                         comparisons = comparisons[[comp]],
                         method = "wilcox.test", label = "p.signif",
                         inherit.aes = TRUE)
    })
}

## 通用：母亲/口腔 → 犊牛生态位 距离数据构建
make_all_distances <- function(st, dist_full) {
  baoji_mothers <- c("Baoji_M0dF", "Baoji_M0dR", "Baoji_M21dF", "Baoji_M21dR")
  gansu_mothers <- c("Gansu_M0dF", "Gansu_M0dR")
  
  dt <- rbindlist(list(
    rbindlist(lapply(baoji_mothers, build_pair_dist, st = st, dist_full = dist_full,
                     g2 = "Baoji_C0dOR", comparison = "Baoji\nMother to Oral",
                     source = "Mother")),
    rbindlist(lapply(baoji_mothers, build_pair_dist, st = st, dist_full = dist_full,
                     g2 = "Baoji_C0dME",
                     comparison = "Baoji\n(Mother or Oral) to Meconium",
                     source = "Mother")),
    build_pair_dist(st, dist_full, "Baoji_C0dOR", "Baoji_C0dME",
                    "Baoji\n(Mother or Oral) to Meconium", source = "Oral"),
    rbindlist(lapply(gansu_mothers, build_pair_dist, st = st, dist_full = dist_full,
                     g2 = "Gansu_C0dME", comparison = "Gansu\nMother to Meconium",
                     source = "Mother"))
  ), fill = TRUE)
  
  dt[, Mother_Group := sub("^[^_]+_", "", Group1)]
  dt[Source == "Oral", Mother_Group := "Oral"]
  dt[, Mother_Group := factor(Mother_Group,
                              levels = c("M21dR", "M21dF", "M0dR", "M0dF", "Oral"))]
  dt[, Comparison := factor(Comparison,
                            levels = c("Baoji\nMother to Oral",
                                       "Baoji\n(Mother or Oral) to Meconium",
                                       "Gansu\nMother to Meconium"))]
  dt
}

DIST_COMPARISONS <- list(
  "Baoji\nMother to Oral" = list(c("M0dF", "M0dR"), c("M21dF", "M21dR"),
                                 c("M0dF", "M21dF"), c("M0dR", "M21dR")),
  "Baoji\n(Mother or Oral) to Meconium" = list(c("M0dF", "M0dR"), c("M21dF", "M21dR"),
                                               c("M0dF", "M21dF"), c("M0dR", "M21dR"),
                                               c("Oral", "M0dF")),
  "Gansu\nMother to Meconium" = list(c("M0dF", "M0dR"))
)

ALPHA_COMPARISONS <- list(
  Baoji = list(c("Baoji_M0dF", "Baoji_C0dME"), c("Baoji_M0dR", "Baoji_C0dME"),
               c("Baoji_M0dF", "Baoji_C0dOR"), c("Baoji_M0dR", "Baoji_C0dOR")),
  Gansu = list(c("Gansu_C0dME", "Gansu_C70dF"), c("Gansu_C70dF", "Gansu_M0dF"),
               c("Gansu_M0dF", "Gansu_M0dR"),  c("Gansu_C70dF", "Gansu_M0dR"),
               c("Gansu_C0dME", "Gansu_M0dF"), c("Gansu_C0dME", "Gansu_M0dR"))
)



# ############################################################################
#  PART I.  病毒 (vOTU)
# ############################################################################

# ----------------------------------------------------------------------------
# Step 1  vOTU 注释表
# ----------------------------------------------------------------------------

viruse_cluster <- fread(file.path(DATA_DIR, "cluster_viruse.tsv"),
                        header = FALSE, sep = "\t",
                        col.names = c("Cluster", "Contig"))
viruse_cluster[, c("Cluster", "Contig") := .(as.character(Cluster), as.character(Contig))]
setkey(viruse_cluster, Cluster)

first_contig <- viruse_cluster[, .(Contig = Contig[1L]), by = Cluster]
setkey(first_contig, Cluster)

anno_v <- fread(file.path(DATA_DIR, "blast_viruse4.tsv"))
anno_v[, c("host_annotation2", "host_Lineage", "lineage") :=
         .(as.character(host_annotation2), as.character(host_Lineage), as.character(lineage))]

log_n <- function(tag, before, after, verb = "removed")
  message(sprintf("  [%s] %d -> %d (%s %d)", tag, before, after, verb, before - after))

# 1.1 排除 Methanobacterium / Methanobrevibacter
n0 <- nrow(anno_v)
anno_v <- anno_v[is.na(host_annotation2) |
                   !(host_annotation2 %in% c("Methanobacterium", "Methanobrevibacter"))]
log_n("Filter 1.1 排除 Methanobacterium/Methanobrevibacter", n0, nrow(anno_v))

# 1.2 排除古菌与真核宿主
n0 <- nrow(anno_v)
bad <- grepl("k__Archaea|k__Eukaryota", anno_v$host_Lineage) |
  grepl("k__Archaea|k__Eukaryota", anno_v$lineage)
anno_v <- anno_v[!bad]
log_n("Filter 1.2 排除 Archaea/Eukaryota", n0, nrow(anno_v))

# 1.3 清除病毒自身注释
idx <- which(grepl("k__Viruse|k__Viridae", anno_v$host_Lineage) |
               grepl("k__Viruse|k__Viridae", anno_v$lineage))
set(anno_v, idx, c("host_annotation2", "host_Lineage"), NA_character_)
message(sprintf("  [Filter 1.3] 清除病毒自身注释: %d 行", length(idx)))

# 1.4 清除 viru / Candidatus
idx <- which(grepl("viru|Candidatus", anno_v$host_annotation2, ignore.case = TRUE))
set(anno_v, idx, "host_annotation2", NA_character_)
message(sprintf("  [Filter 1.4] 清除 viru/Candidatus: %d 行", length(idx)))

# 1.5 CheckV 质量筛选
checkv <- rbindlist(lapply(file.path(DATA_DIR, c("quality_summary_checkv.tsv",
                                                 "quality_summary_checkv2.tsv")), fread),
                    fill = TRUE)
checkv <- unique(checkv, by = "contig_id")

checkv_all <- copy(checkv)
checkv_all[, Passed := contamination < MAX_CONTAM & checkv_quality %in% QUALITY_KEEP]
checkv_all[is.na(checkv_quality) | checkv_quality == "", checkv_quality := "Not-filtered"]
log_n("Filter 1.5 CheckV", nrow(checkv),
      checkv[contamination < MAX_CONTAM & checkv_quality %in% QUALITY_KEEP, .N])
checkv <- checkv[contamination < MAX_CONTAM & checkv_quality %in% QUALITY_KEEP,
                 .(contig_id, checkv_quality)]
setkey(checkv, contig_id)

# 1.6 仅保留通过 CheckV 的 contig
n0 <- nrow(anno_v)
anno_v <- anno_v[Contig_name %in% checkv$contig_id]
log_n("Filter 1.6 CheckV 过滤", n0, nrow(anno_v))
anno_v <- unique(anno_v, by = "Contig_name")
setkey(anno_v, Contig_name)

# 1.7 合并 聚类 + 注释 + 质量
viruse_anno_full <- merge(viruse_cluster, anno_v, by.x = "Cluster", by.y = "Contig_name")
viruse_anno_full <- merge(viruse_anno_full, checkv,
                          by.x = "Cluster", by.y = "contig_id", all.x = TRUE)
message(sprintf("  [Filter 1.7] 合并后: %d vOTUs", nrow(viruse_anno_full)))

# ============================================================================
# extract_tax_level: 从 lineage 字符串中提取指定分类水平
#   支持 "k__X;p__Y;..." 带前缀格式 与 无前缀按位置格式
#   支持 ; 或 | 分隔
# ============================================================================
extract_tax_level <- function(lineage, level) {
  
  prefix_map <- c(kingdom = "k", superkingdom = "k", domain = "d",
                  phylum = "p", class = "c", order = "o",
                  family = "f", genus = "g", species = "s",
                  strain = "t")
  key <- tolower(as.character(level))[1]
  pfx <- prefix_map[[key]]
  if (is.null(pfx)) stop("不支持的分类水平: ", level)
  
  s <- as.character(lineage)
  
  ## ---- (1) 带前缀匹配：p__Firmicutes ----
  pat <- paste0("(?:^|[;|]\\s*)", pfx, "__\\s*([^;|]*)")
  m   <- regmatches(s, regexec(pat, s, ignore.case = TRUE, perl = TRUE))
  val <- vapply(m, function(x) if (length(x) >= 2L) trimws(x[2L]) else NA_character_,
                character(1))
  
  ## ---- (2) 无前缀时按位置取（k=1, p=2, c=3, o=4, f=5, g=6, s=7）----
  need_fb <- is.na(val) & !is.na(s) & nzchar(s)
  if (any(need_fb)) {
    pos <- which(c("k", "p", "c", "o", "f", "g", "s") == pfx)
    if (length(pos)) {
      fb <- vapply(s[need_fb], function(one) {
        parts <- trimws(unlist(strsplit(one, "[;|]", perl = TRUE)))
        parts <- parts[nzchar(parts)]
        if (length(parts) >= pos) sub("^[a-z]__\\s*", "", parts[pos], ignore.case = TRUE)
        else NA_character_
      }, character(1), USE.NAMES = FALSE)
      val[need_fb] <- fb
    }
  }
  
  ## ---- (3) 清洗占位值 → NA ----
  bad <- c("", "-", "--", "na", "n/a", "none", "null", "unknown", "uncultured",
           "unclassified", "metagenome", "environmental", "environmental samples",
           "unidentified", "no rank", "incertae sedis", "0")
  val[is.na(val)] <- ""
  val[tolower(val) %in% bad] <- NA_character_
  val[grepl("^(unclassified|uncultured|metagenom)", val, ignore.case = TRUE)] <- NA_character_
  val[!nzchar(val)] <- NA_character_
  
  val
}
# 1.8 解析宿主 lineage
for (lv in c("Kingdom", "Phylum", "Class", "Order", "Family")) {
  set(viruse_anno_full, NULL, paste0("host_", lv),
      extract_tax_level(viruse_anno_full$host_Lineage, lv))
}
set(viruse_anno_full, NULL, "host_Genus", viruse_anno_full$host_annotation2)

# 1.9 标准注释表
viruse_anno_clean <- unique(viruse_anno_full[, .(
  OTU = Cluster, Kingdom = "k__Viruses",
  Phylum = NA_character_, Class = NA_character_, Order = NA_character_,
  Family = annotation_pre, Genus = NA_character_, Species = NA_character_,
  checkv_quality, host_Kingdom, host_Phylum, host_Class, host_Order,
  host_Family, host_Genus)], by = "OTU")

# 1.10 长度
viruse_anno_clean <- merge(viruse_anno_clean, first_contig,
                           by.x = "OTU", by.y = "Cluster", all.x = TRUE, sort = FALSE)
set(viruse_anno_clean, NULL, "length", extract_contig_length(viruse_anno_clean$Contig))
viruse_anno_clean[is.na(length), length := extract_contig_length(OTU[is.na(length)])]
viruse_anno_clean[, Contig := NULL]

n0 <- nrow(viruse_anno_clean)
viruse_anno_clean <- viruse_anno_clean[is.finite(length) & length > 5000]
log_n("Filter 1.8 length > 0", n0, nrow(viruse_anno_clean))

# 1.11 Family 缺失
viruse_anno_clean[is.na(Family) | Family == "", Family := "Unclassified"]
message(sprintf("  [Filter 1.9] Unclassified Family: %d vOTUs",
                viruse_anno_clean[Family == "Unclassified", .N]))
setkey(viruse_anno_clean, OTU)

save_tsv(viruse_anno_clean, "clean_votu_annotation.tsv")

message(sprintf("  \u2713 Final clean vOTU annotation: %d vOTUs", nrow(viruse_anno_clean)))



# ----------------------------------------------------------------------------
# Step 1.13  CheckV 质量分布图（完整度 高 → 低，颜色 绿 → 黄）
# ----------------------------------------------------------------------------

## 1) 质量等级顺序：按完整度从高到低
QUALITY_ORDER <- c("Complete", "High-quality", "Medium-quality",
                   "Low-quality", "Not-filtered")

## 2) 绿 → 黄 连续配色（Not-filtered 固定灰色）
build_quality_palette <- function(levs) {
  real <- levs[levs != "Not-filtered"]
  pal  <- colorRampPalette(c("#145A32",   # 深绿  Complete
                             "#2E8B57",   # 海绿  High
                             "#8FBF5A",   # 黄绿  Medium
                             "#D9D24A",   # 浅黄
                             "#F2C14E"))(length(real))   # 暖黄  Low
  names(pal) <- real
  if ("Not-filtered" %in% levs) pal <- c(pal, "Not-filtered" = "#BDBDBD")
  pal[levs]
}

## 3) 读入 + 去重 + 标注是否通过
checkv <- rbindlist(lapply(file.path(DATA_DIR,
                                     c("quality_summary_checkv.tsv",
                                       "quality_summary_checkv2.tsv")), fread),
                    fill = TRUE)
checkv <- unique(checkv, by = "contig_id")

checkv_all <- copy(checkv)
checkv_all[, Passed := contamination < MAX_CONTAM & checkv_quality %in% QUALITY_KEEP]
checkv_all[is.na(checkv_quality) | checkv_quality == "", checkv_quality := "Not-filtered"]

## 4) 只保留最终 vOTU，再因子化（顺序 = 完整度高→低）
checkv_all <- checkv_all[contig_id %in% viruse_anno_clean$OTU]
checkv_all[, checkv_quality := factor(checkv_quality,
                                      levels = intersect(QUALITY_ORDER,
                                                         unique(checkv_quality)))]
checkv_all <- checkv_all[!is.na(checkv_quality)]
checkv_all[, checkv_quality := droplevels(checkv_quality)]

## 5) 汇总
qual_sum <- checkv_all[, .(N                  = .N,
                           Mean_Completeness  = round(mean(completeness,   na.rm = TRUE), 1),
                           Mean_Contamination = round(mean(contamination,  na.rm = TRUE), 2)),
                       by = checkv_quality]
qual_sum[, checkv_quality := factor(checkv_quality,
                                    levels = levels(checkv_all$checkv_quality))]
setorder(qual_sum, checkv_quality)                       # 按完整度高→低
qual_sum[, Pct   := round(N / sum(N) * 100, 1)]
qual_sum[, Label := sprintf("%d\n(%.1f%%)", N, Pct)]

# X 轴标签：等级名 + 平均完整度
qual_sum[, X_Label := ifelse(is.na(Mean_Completeness),
                             as.character(checkv_quality),
                             sprintf("%s\n(comp %s%%)", checkv_quality, Mean_Completeness))]
x_lab <- setNames(qual_sum$X_Label, as.character(qual_sum$checkv_quality))

print(qual_sum)

## 6) 绘图
QUALITY_COLOR <- build_quality_palette(levels(qual_sum$checkv_quality))

p_checkv <- ggplot(qual_sum, aes(x = checkv_quality, y = N, fill = checkv_quality)) +
  geom_col(width = 0.68, color = "grey25", linewidth = 0.3) +
  geom_text(aes(label = Label), vjust = -0.25, size = 3,
            lineheight = 0.85, color = "grey15") +
  scale_x_discrete(limits = levels(qual_sum$checkv_quality),
                   labels = x_lab, drop = FALSE) +
  scale_fill_manual(values = QUALITY_COLOR, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.20))) +
  # labs(x = NULL, y = "Number of vOTUs",
  #      title = "CheckV quality of retained vOTUs",
  #      subtitle = sprintf("n = %d  |  ordered by genome completeness (high \u2192 low)",
  #                         nrow(checkv_all))) +
  theme_bw(base_size = 9) +
  theme(legend.position    = "none",
        plot.title         = element_text(size = 10, face = "bold", hjust = 0.5),
        plot.subtitle      = element_text(size = 7.5, hjust = 0.5, color = "grey35"),
        axis.text.x        = element_text(size = 7.5, lineheight = 0.9),
        axis.text.y        = element_text(size = 8),
        axis.title.y       = element_text(size = 9),
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        plot.margin        = margin(8, 10, 8, 8, unit = "pt"))

save_pdf(p_checkv, "Fig1_CheckV_Quality.pdf", 8, 6)
save_tsv(qual_sum, "CheckV_quality_distribution.tsv")
message(sprintf("  \u2713 CheckV distribution: %d vOTUs across %d quality tiers",
                nrow(checkv_all), nrow(qual_sum)))


anno_v2 <- anno_v[Contig_name %in% viruse_anno_clean$OTU]
anno_v2 <- merge(anno_v2,checkv_all[,c("contig_id","checkv_quality")],by.x="Contig_name",by.y="contig_id")
head(anno_v2)
write.table(anno_v2, "../data/blast_viruse.tsv",sep = "\t",quote = F,row.names = F)



library(data.table)
library(stringr)
anno_v2[, Lineage_Genus := str_extract(lineage, "(?<=g__)[^;]+")]
dt_valid <- anno_v2[!is.na(host_annotation) & host_annotation != "" & 
                      !is.na(Lineage_Genus) & Lineage_Genus != ""]
dt_valid <- dt_valid[-grep("unclassified",dt_valid$Lineage_Genus),]
dt_valid <- dt_valid[-grep("g__unclassified",dt_valid$host_Lineage),]
dt_valid <- dt_valid[is.na(dt_valid$host_Lineage)==F,]
n_consistent <- sum(dt_valid$host_annotation == dt_valid$Lineage_Genus)
n_inconsistent <- nrow(dt_valid) - n_consistent

# ----------------------------------------------------------------------------
# Step 2  vOTU RPKM 矩阵
# ----------------------------------------------------------------------------


# 2.1 元数据
meta <- fread(file.path(DATA_DIR, "meta.tsv"))
n_meta_total <- nrow(meta)
meta <- meta[!Group %in% EXCLUDE_GROUP]
message(sprintf("  [Filter 2.1] Meta: %d -> %d (排除 %s)",
                n_meta_total, nrow(meta), paste(EXCLUDE_GROUP, collapse = "/")))

id_col <- intersect(c("Sample_ID", "sample_id", "SampleID", "sample", "Sample", "ID"),
                    names(meta))[1]
if (is.na(id_col)) stop("meta.tsv 中找不到样本 ID 列")
meta[[id_col]] <- trimws(as.character(meta[[id_col]]))
setnames(meta, id_col, "Sample_ID", skip_absent = TRUE)
valid_samples <- meta$Sample_ID

# 2.2 counts
counts_dt <- fread(file.path(DATA_DIR, "fig1/merged_abundance_matrix.tsv"))
setnames(counts_dt, 1L, "OTU")
counts_dt[, OTU := as.character(OTU)]
all_sample_cols <- setdiff(names(counts_dt), "OTU")
samp_cols <- intersect(valid_samples, all_sample_cols)
n_otu_total <- nrow(counts_dt)
message(sprintf("  [Filter 2.2] 丰度矩阵: %d vOTUs x %d samples", n_otu_total, length(all_sample_cols)))

# 2.3 仅保留有注释的 OTU
counts_dt <- counts_dt[OTU %in% viruse_anno_clean$OTU, c("OTU", samp_cols), with = FALSE]
log_n("Filter 2.3 仅保留注释 vOTU", n_otu_total, nrow(counts_dt))

# 2.4 长度筛选
anno <- viruse_anno_clean[length >= MIN_CONTIG_LEN]
n0 <- nrow(counts_dt)
counts_dt <- counts_dt[OTU %in% anno$OTU]
log_n(sprintf("Filter 2.4 length >= %d bp", MIN_CONTIG_LEN), n0, nrow(counts_dt))

# 2.5 RPKM
cnt_mat <- as.matrix(counts_dt[, samp_cols, with = FALSE])
rownames(cnt_mat) <- counts_dt$OTU
storage.mode(cnt_mat) <- "double"
rm(counts_dt); invisible(gc())

anno    <- anno[match(rownames(cnt_mat), OTU)]
len_kb  <- anno$length / 1000
depth_m <- meta[match(colnames(cnt_mat), Sample_ID), num_seqs] / 1e6
ok_r <- is.finite(len_kb)  & len_kb  > 0
ok_c <- is.finite(depth_m) & depth_m > 0

rpkm <- matrix(0, nrow(cnt_mat), ncol(cnt_mat), dimnames = dimnames(cnt_mat))
if (any(ok_r) && any(ok_c)) {
  rpkm[ok_r, ok_c] <- sweep(sweep(cnt_mat[ok_r, ok_c, drop = FALSE],
                                  1, len_kb[ok_r], "/"), 2, depth_m[ok_c], "/")
}
rm(cnt_mat); invisible(gc())

# 2.6 RPKM 阈值
n0 <- nrow(rpkm)
keep <- row_max(rpkm) >= MIN_RPKM
rpkm <- rpkm[keep, , drop = FALSE]; anno <- anno[keep]
log_n(sprintf("Filter 2.5 max RPKM >= %g", MIN_RPKM), n0, nrow(rpkm))

# 2.7 全 0 行
n0 <- nrow(rpkm)
keep <- rowSums(rpkm>0) > 0
rpkm <- rpkm[keep, , drop = FALSE]; anno <- anno[keep]
log_n("Filter 2.6 去除全 0 行", n0, nrow(rpkm))

anno_W <- anno_v[anno_v$Contig_name%in%anno$OTU,]%>%
  select(Contig_name, Sample, host_annotation2,annotation_pre, predicted_label, Individuality, Time, Type, 
         Paired_number, Age, Group, Area)
save_tsv(data.table(OTU = rownames(rpkm), as.data.table(rpkm)),
         "clean_votu_rpkm_filtered.tsv")
message(sprintf("  \u2713 Final: %d vOTUs x %d samples", nrow(rpkm), ncol(rpkm)))





counts_dt <- fread(file.path(DATA_DIR, "fig1/merged_abundance_matrix.tsv"))
setnames(counts_dt, 1L, "OTU")
counts_dt[, OTU := as.character(OTU)]
all_sample_cols <- setdiff(names(counts_dt), "OTU")
samp_cols <- intersect(valid_samples, all_sample_cols)
n_otu_total <- nrow(counts_dt)
message(sprintf("  [Filter 2.2] 丰度矩阵: %d vOTUs x %d samples", n_otu_total, length(all_sample_cols)))



# 2.5 RPKM
cnt_mat <- as.matrix(counts_dt[, samp_cols, with = FALSE])
rownames(cnt_mat) <- counts_dt$OTU
storage.mode(cnt_mat) <- "double"
rm(counts_dt); invisible(gc())

anno    <- anno[match(rownames(cnt_mat), OTU)]
len_kb  <- anno$length / 1000
depth_m <- meta[match(colnames(cnt_mat), Sample_ID), num_seqs] / 1e6
ok_r <- is.finite(len_kb)  & len_kb  > 0
ok_c <- is.finite(depth_m) & depth_m > 0

rpkm <- matrix(0, nrow(cnt_mat), ncol(cnt_mat), dimnames = dimnames(cnt_mat))
if (any(ok_r) && any(ok_c)) {
  rpkm[ok_r, ok_c] <- sweep(sweep(cnt_mat[ok_r, ok_c, drop = FALSE],
                                  1, len_kb[ok_r], "/"), 2, depth_m[ok_c], "/")
}
rm(cnt_mat); invisible(gc())

# ----------------------------------------------------------------------------
# Step 3  sample_table / tax_table / microeco 对象
#         （原 Step 7.1 的因子标准化已合并到此处）
# ----------------------------------------------------------------------------


sample_table <- as.data.table(meta)[Sample_ID %in% colnames(rpkm)]
sample_table <- sample_table[match(colnames(rpkm), Sample_ID)]
setDT(sample_table)

# Group_Area
if (!"Group_Area" %in% names(sample_table)) {
  sample_table[, Group_Area := paste0(Area, "_", Group)]
}

# Time
sample_table[, Time := fcase(Time %in% c("0", "d0"), "d0",
                             Time %in% c("-21", "21"), "d-21",
                             default = "d70")]
sample_table[, Time := factor(Time, levels = c("d-21", "d0", "d70"))]

# Niche
sample_table[, Niche := fcase(Type %in% c("Oral"),              "Oral",
                              Type %in% c("fecal", "FE"),       "Fecal",
                              Type %in% c("rumen", "RF"),       "Rumen",
                              default = as.character(Type))]
sample_table[, Niche := factor(Niche, levels = c("Oral", "Fecal", "Rumen"))]

# Age / Area / Group_Area
sample_table[, Age  := factor(Age,  levels = c("Calf", "Mother", "Adult"))]
sample_table[, Area := factor(Area, levels = c("Baoji", "Gansu"))]
sample_table[, Group_Area := factor(Group_Area,
                                    levels = intersect(GROUP_LEVELS, unique(Group_Area)))]
rownames(sample_table) <- sample_table$Sample_ID

# tax_table（行序与 rpkm 一致）
tax_table <- as.data.frame(anno[match(rownames(rpkm), OTU), TAX_LEVELS, with = FALSE])
rownames(tax_table) <- rownames(rpkm)

message(sprintf("  \u2713 sample_table: %d samples | tax_table: %d OTUs x %d ranks",
                nrow(sample_table), nrow(tax_table), ncol(tax_table)))


# ----------------------------------------------------------------------------
# Step 4  病毒科水平丰度（排除 Unclassified）
# ----------------------------------------------------------------------------
# ============================================================================
# make_microeco: 构建 microeco 兼容对象（强制 data.frame + 对齐 + 去 NA）
# ============================================================================
make_microeco <- function(rpkm_mat, anno_dt, sample_table_dt,
                          keep_idx = NULL, tax_levels = TAX_LEVELS) {
  if (is.null(keep_idx)) keep_idx <- seq_len(nrow(rpkm_mat))
  
  otu_ids <- rownames(rpkm_mat)[keep_idx]
  samp_ids <- colnames(rpkm_mat)
  
  ## 1) otu_table：matrix / data.frame，行 = OTU，列 = 样本
  otu_df <- as.data.frame(rpkm_mat[keep_idx, , drop = FALSE])
  rownames(otu_df) <- otu_ids
  
  ## 2) tax_table：行名必须与 otu_table 完全一致，且不能有 NA
  tax_df <- as.data.frame(anno_dt[match(otu_ids, anno_dt$OTU), tax_levels, with = FALSE])
  rownames(tax_df) <- otu_ids
  tax_df <- tax_df[!is.na(match(rownames(tax_df), otu_ids)), , drop = FALSE]
  for (cc in names(tax_df)) {
    v <- as.character(tax_df[[cc]])
    v[is.na(v) | !nzchar(v)] <- "Unclassified"
    tax_df[[cc]] <- v
  }
  
  ## 3) sample_table：必须是 data.frame，行名 = Sample_ID
  samp_df <- as.data.frame(sample_table_dt)
  rownames(samp_df) <- samp_df$Sample_ID
  samp_df <- samp_df[intersect(samp_ids, samp_df$Sample_ID), , drop = FALSE]
  for (cc in intersect(c("Group", "Area", "Niche", "Age", "Time", "Group_Area"),
                       names(samp_df))) {
    samp_df[[cc]] <- as.character(samp_df[[cc]])
  }
  
  ## 4) 关键：otu_table 列名与 sample_table 行名严格对齐
  shared <- intersect(colnames(otu_df), rownames(samp_df))
  otu_df <- otu_df[, shared, drop = FALSE]
  samp_df <- samp_df[shared, , drop = FALSE]
  otu_df <- otu_df[rowSums(otu_df) > 0, , drop = FALSE]
  
  ## 5) 三个对象都不能是 data.table
  stopifnot(!inherits(otu_df, "data.table"),
            !inherits(tax_df, "data.table"),
            !inherits(samp_df, "data.table"))
  
  microtable$new(sample_table = samp_df,
                 otu_table    = otu_df,
                 tax_table    = tax_df)
}


keep_idx <- which(anno$Family != "Unclassified")
message(sprintf("  [Filter 4.1] 排除 Unclassified: %d -> %d vOTUs",
                nrow(anno), length(keep_idx)))

mt_plot <- make_microeco(rpkm, anno, sample_table, keep_idx = keep_idx)

message(sprintf("  microeco 对象: %d OTUs x %d samples | tax ranks = %s",
                nrow(mt_plot$otu_table), ncol(mt_plot$otu_table),
                paste(colnames(mt_plot$tax_table), collapse = "/")))

t1 <- trans_abund$new(dataset = mt_plot, taxrank = "Family",
                      ntaxa = 5, groupmean = "Group")

p_abund <- t1$plot_bar(others_color = "grey70", xtext_keep = TRUE,
                       legend_text_italic = FALSE,
                       color_values = BAR_COLORS[seq_len(5)]) +
  theme(axis.title.y    = element_text(size = 8),
        axis.text.x     = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
        axis.text.y     = element_text(size = 7),
        legend.text     = element_text(size = 7),
        legend.title    = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"))

save_pdf(p_abund, "Fig1_Viral_Abundance.pdf", 8, 6)


# ----------------------------------------------------------------------------
# Step 5  病毒科 ↔ 宿主门 连接图
# ----------------------------------------------------------------------------


vh <- unique(anno[!is.na(Family) & !is.na(host_Phylum) & host_Phylum != "",
                  .(OTU, virus_family = Family, host_phylum = host_Phylum)])
message(sprintf("  [Filter 5.1] 同时具备 Family 与 host_Phylum: %d vOTUs (from %d)",
                uniqueN(vh$OTU), nrow(anno)))

host_rank  <- vh[, .(n_otu = uniqueN(OTU)), by = host_phylum][order(-n_otu)]
top_host   <- head(host_rank$host_phylum, 5)
host_order <- c(top_host)
vh[, host_phylum_plot := fifelse(host_phylum %in% top_host, host_phylum, "Others")]

fam_rank <- vh[, .(total_otu = uniqueN(OTU)), by = virus_family][order(-total_otu)]
target <- c("Siphoviridae", "Myoviridae", "Podoviridae", "Steigviridae", "Herelleviridae")
avail  <- intersect(target, fam_rank$virus_family)
virus_order <- if (length(avail) >= 5) avail else head(unique(c(avail, fam_rank$virus_family)), 5)

link <- vh[virus_family %in% virus_order,
           .(n_otu = uniqueN(OTU)), by = .(virus_family, host_phylum_plot)]
link[, y_left  := length(virus_order) + 1 - match(virus_family, virus_order)]
link[, y_right := length(host_order)  + 1 - match(host_phylum_plot, host_order)]

left_nodes  <- unique(link[, .(label = virus_family,     y = y_left)])
right_nodes <- unique(link[, .(label = host_phylum_plot, y = y_right)])
host_palette <- setNames(c(brewer.pal(8, "Set2")[seq_along(top_host)], "#BFBFBF"), host_order)

p_host <- ggplot(link) +
  geom_segment(aes(x = 1, xend = 2, y = y_left, yend = y_right, color = host_phylum_plot),
               linewidth = 0.8, alpha = 0.7, lineend = "round") +
  geom_point(data = left_nodes, aes(x = 1, y = y),
             size = 3.5, shape = 21, fill = "white", color = "black", stroke = 0.7) +
  geom_point(data = right_nodes, aes(x = 2, y = y, color = label),
             size = 3.5, shape = 21, fill = "white", stroke = 0.7) +
  geom_text(data = left_nodes, aes(x = 1, y = y, label = label),
            hjust = 1.15, size = 3, fontface = "italic") +
  geom_text(data = right_nodes, aes(x = 2, y = y, label = label),
            hjust = -0.15, size = 3, fontface = "italic") +
  scale_color_manual(values = host_palette, name = "Host Phylum") +
  scale_x_continuous(limits = c(0.55, 2.50), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.5, max(length(virus_order), length(host_order)) + 0.5),
                     expand = c(0, 0)) +
  theme_void(base_size = 14) +
  theme(legend.position = "none",
        plot.margin = margin(10, 120, 10, 120, unit = "pt"))

save_pdf(p_host, "Fig1_Virus_Host_Phylum.pdf", 14, 10)


# ----------------------------------------------------------------------------
# Step 6  病毒 Alpha 多样性
# ----------------------------------------------------------------------------


alpha_data <- data.table(Sample_ID = colnames(rpkm),
                         shannon   = vegan::diversity(t(rpkm), index = "shannon"))
alpha_data <- merge(alpha_data,
                    as.data.table(sample_table)[, .(Sample_ID, Group_Area, Area, Group)],
                    by = "Sample_ID")
set(alpha_data, NULL, "Group_Area", factor(alpha_data$Group_Area, levels = GROUP_LEVELS))

p_alpha <- plot_alpha_box(alpha_data, ALPHA_COMPARISONS,
                          group_color[levels(alpha_data$Group_Area)])
save_pdf(p_alpha, "Fig1_Shannon_Diversity.pdf", 14, 5)
save_tsv(alpha_data, "Viral_alpha_diversity.tsv")
message("  \u2713 Viral Shannon diversity plot saved")


# ----------------------------------------------------------------------------
# Step 7  病毒 Beta 多样性
# ----------------------------------------------------------------------------


beta_v <- run_beta(rpkm, sample_table, tag = "Viral ", w = 10, h = 8)
save_pdf(beta_v$plot, "Fig1_PCoA.pdf", 10, 8)
save_tsv(beta_v$sig,  "Beta_diversity_significance.tsv")
message(sprintf("  \u2713 PERMANOVA R2 = %.3f, P %s | Betadisper P = %.3f",
                beta_v$r2_permanova, fmt_p(beta_v$p_permanova), beta_v$p_disp))


# ----------------------------------------------------------------------------
# Step 7.2  病毒方差分解
# ----------------------------------------------------------------------------


vp_v <- run_variance_partition(beta_v$dist_full, beta_v$st, label = "Viral")
save_tsv(vp_v$single, "Variance_single_factor.tsv")
save_tsv(vp_v$multi,  "Variance_multi_factor.tsv")
save_tsv(vp_v$cohort, "Variance_cohort_specific.tsv")

# ----------------------------------------------------------------------------
# Step 7.3-7.5  病毒组间距离
# ----------------------------------------------------------------------------

all_distances <- make_all_distances(beta_v$st, beta_v$dist_full)

dist_summary <- all_distances[, .(N = .N,
                                  Mean_Distance   = round(mean(Distance),   4),
                                  SD_Distance     = round(sd(Distance),     4),
                                  Median_Distance = round(median(Distance), 4)),
                              by = .(Comparison, Mother_Group)]
print(dist_summary)
save_tsv(dist_summary,  "Inter_group_distances_summary.tsv")
save_tsv(all_distances, "Inter_group_distances_all.tsv")

p_dist_viruse <- plot_dist_box(all_distances, DIST_COMPARISONS,
                        "Viral Community Distance: Maternal Groups to Calf Niches")
save_pdf(p_dist_viruse, "Fig1_Inter_Group_Distances.pdf", 8, 4)

dist_tests <- rbindlist(lapply(split(all_distances, all_distances$Comparison), function(d) {
  if (d[, uniqueN(Mother_Group)] < 2) return(NULL)
  kw  <- kruskal.test(Distance ~ Mother_Group, data = d)
  out <- list(data.table(Test = "Kruskal-Wallis",
                         Statistic = round(unname(kw$statistic), 3),
                         P_value = kw$p.value))
  if (kw$p.value < 0.05) {
    pw <- pairwise.wilcox.test(d$Distance, d$Mother_Group, p.adjust.method = "fdr")$p.value
    ij <- which(!is.na(pw), arr.ind = TRUE)
    out[[2]] <- data.table(Test = sprintf("Wilcoxon: %s vs %s",
                                          rownames(pw)[ij[, 1]], colnames(pw)[ij[, 2]]),
                           Statistic = NA_real_, P_value = pw[ij])
  }
  rbindlist(out, fill = TRUE)
}), idcol = "Comparison", fill = TRUE)
print(dist_tests)
save_tsv(dist_tests, "Inter_group_distance_tests.tsv")


# ############################################################################
#  PART II.  细菌 (bOTU)
# ############################################################################

# ----------------------------------------------------------------------------
# Step B1  细菌数据预处理与 RPKM
# ----------------------------------------------------------------------------


# B1.1 聚类表 + 长度 > MIN_CONTIG_LEN
cluster_bact <- fread(file.path(DATA_DIR, "cluster_bact.tsv"),
                      header = FALSE, sep = "\t",
                      col.names = c("Cluster", "Contig"))
cluster_bact[, c("Cluster", "Contig") := .(as.character(Cluster), as.character(Contig))]
setkey(cluster_bact, Cluster)

first_contig_bact <- cluster_bact[, .(Contig = Contig[1L]), by = Cluster]
setkey(first_contig_bact, Cluster)

cluster_bact <- cluster_bact[extract_contig_length(Cluster) > MIN_CONTIG_LEN]
message(sprintf("  [Filter B1.1] length > %d bp clusters: %d",
                MIN_CONTIG_LEN, uniqueN(cluster_bact$Cluster)))

# B1.2 丰度矩阵（自动识别 ID 列，切勿先删列）
bact_abundance <- fread(file.path(DATA_DIR, "merged_abundance_matrix_bact222.tsv"))
id_col_b <- intersect(c("contig", "OTU", "Cluster", "cluster"), names(bact_abundance))[1]
if (is.na(id_col_b)) id_col_b <- names(bact_abundance)[1]
setnames(bact_abundance, id_col_b, "OTU", skip_absent = TRUE)
bact_abundance[, OTU := as.character(OTU)]

n_bact_total <- nrow(bact_abundance)
bact_abundance <- bact_abundance[OTU %in% cluster_bact$Cluster]
log_n("Filter B1.2 按 cluster 过滤", n_bact_total, nrow(bact_abundance))

# B1.3 注释表
blast_bact <- fread(file.path(DATA_DIR, "blast_bact.tsv"))
blast_bact <- unique(blast_bact[Contig_name %in% cluster_bact$Cluster], by = "Contig_name")
setkey(blast_bact, Contig_name)
message(sprintf("  [Filter B1.3] 细菌注释: %d", nrow(blast_bact)))

# B1.4 组装标准注释表
pick_tax_col <- function(dt, lv) {
  cands <- c(paste0("tax_", lv), paste0("tax_", tolower(lv)),
             lv, tolower(lv), paste0("taxon_", lv), paste0("Tax_", lv))
  hit <- intersect(cands, names(dt))
  if (length(hit)) hit[1] else NA_character_
}

bact_anno_full <- merge(cluster_bact, blast_bact, by.x = "Cluster", by.y = "Contig_name")
bact_anno_full <- unique(bact_anno_full, by = "Cluster")

mi <- match(bact_anno_full$Cluster, blast_bact$Contig_name)
bact_anno_clean <- data.table(OTU = bact_anno_full$Cluster)
for (lv in TAX_LEVELS) {
  cc <- pick_tax_col(blast_bact, lv)
  set(bact_anno_clean, NULL, lv,
      if (is.na(cc)) NA_character_ else as.character(blast_bact[[cc]][mi]))
}
bact_anno_clean <- unique(bact_anno_clean, by = "OTU")

# B1.5 长度
bact_anno_clean <- merge(bact_anno_clean, first_contig_bact,
                         by.x = "OTU", by.y = "Cluster", all.x = TRUE, sort = FALSE)
set(bact_anno_clean, NULL, "length", extract_contig_length(bact_anno_clean$Contig))
bact_anno_clean[is.na(length), length := extract_contig_length(OTU[is.na(length)])]
bact_anno_clean[, Contig := NULL]

n0 <- nrow(bact_anno_clean)
bact_anno_clean <- bact_anno_clean[is.finite(length) & length > 5000]
log_n("Filter B1.4 length > 0", n0, nrow(bact_anno_clean))

bact_anno_clean[is.na(Family) | Family == "", Family := "Unclassified"]
message(sprintf("  [Filter B1.5] Unclassified Family: %d bOTUs",
                bact_anno_clean[Family == "Unclassified", .N]))
setkey(bact_anno_clean, OTU)

# B1.6 样本列 + counts
samp_cols_bact <- intersect(meta$Sample_ID, setdiff(names(bact_abundance), "OTU"))
bact_counts <- bact_abundance[#OTU %in% bact_anno_clean$OTU,
                              ,c("OTU", samp_cols_bact), with = FALSE]
message(sprintf("  [Filter B1.6] counts: %d bOTUs x %d samples",
                nrow(bact_counts), length(samp_cols_bact)))

# B1.7 RPKM
cnt_mat_bact <- as.matrix(bact_counts[, samp_cols_bact, with = FALSE])
rownames(cnt_mat_bact) <- bact_counts$OTU
storage.mode(cnt_mat_bact) <- "double"
rm(bact_counts); invisible(gc())

anno_bact    <- bact_anno_clean[match(rownames(cnt_mat_bact), OTU)]
len_kb_bact  <- anno_bact$length / 1000
depth_m_bact <- meta[match(colnames(cnt_mat_bact), Sample_ID), num_seqs] / 1e6
ok_r_b <- is.finite(len_kb_bact)  & len_kb_bact  > 0
ok_c_b <- is.finite(depth_m_bact) & depth_m_bact > 0

rpkm_bact <- matrix(0, nrow(cnt_mat_bact), ncol(cnt_mat_bact), dimnames = dimnames(cnt_mat_bact))
if (any(ok_r_b) && any(ok_c_b)) {
  rpkm_bact[ok_r_b, ok_c_b] <- sweep(sweep(cnt_mat_bact[ok_r_b, ok_c_b, drop = FALSE],
                                           1, len_kb_bact[ok_r_b], "/"),
                                     2, depth_m_bact[ok_c_b], "/")
}
rm(cnt_mat_bact); invisible(gc())

# B1.8 过滤
n0 <- nrow(rpkm_bact)
keep_bact <- row_max(rpkm_bact) >= MIN_RPKM
rpkm_bact <- rpkm_bact[keep_bact, , drop = FALSE]; anno_bact <- anno_bact[keep_bact]
log_n(sprintf("Filter B1.7 max RPKM >= %g", MIN_RPKM), n0, nrow(rpkm_bact))

n0 <- nrow(rpkm_bact)
keep_bact <- rowSums(rpkm_bact) > 0
rpkm_bact <- rpkm_bact[keep_bact, , drop = FALSE]; anno_bact <- anno_bact[keep_bact]
log_n("Filter B1.8 去除全 0 行", n0, nrow(rpkm_bact))

save_tsv(anno_bact, "clean_bact_annotation_filtered.tsv")
save_tsv(data.table(OTU = rownames(rpkm_bact), as.data.table(rpkm_bact)),
         "clean_bact_rpkm_filtered.tsv")
message(sprintf("  \u2713 Final bacterial matrix: %d OTUs x %d samples",
                nrow(rpkm_bact), ncol(rpkm_bact)))


# ----------------------------------------------------------------------------
# Step B2  细菌 Alpha 多样性
# ----------------------------------------------------------------------------


alpha_data_bact <- data.table(Sample_ID = colnames(rpkm_bact),
                              shannon   = vegan::diversity(t(rpkm_bact), index = "shannon"))
alpha_data_bact <- merge(alpha_data_bact,
                         as.data.table(sample_table)[, .(Sample_ID, Group_Area, Area, Group)],
                         by = "Sample_ID")
set(alpha_data_bact, NULL, "Group_Area",
    factor(alpha_data_bact$Group_Area, levels = GROUP_LEVELS))

p_alpha_bact <- plot_alpha_box(alpha_data_bact, ALPHA_COMPARISONS,
                               group_color[levels(alpha_data_bact$Group_Area)])
save_pdf(p_alpha_bact, "Fig1_Bact_Shannon_Diversity.pdf", 14, 5)
save_tsv(alpha_data_bact, "Bact_alpha_diversity.tsv")
message("  \u2713 Bacterial Shannon diversity plot saved")

# ----------------------------------------------------------------------------
# Step B3  细菌 Beta 多样性
# ----------------------------------------------------------------------------


beta_b <- run_beta(rpkm_bact, sample_table, tag = "Bacterial ", w = 10, h = 8)
save_pdf(beta_b$plot, "Fig1_Bact_PCoA.pdf", 10, 8)
save_tsv(beta_b$sig,  "Bact_Beta_diversity_significance.tsv")
message(sprintf("  \u2713 PERMANOVA R2 = %.3f, P %s | Betadisper P = %.3f",
                beta_b$r2_permanova, fmt_p(beta_b$p_permanova), beta_b$p_disp))


# ----------------------------------------------------------------------------
# Step B4  细菌方差分解
# ----------------------------------------------------------------------------

vp_b <- run_variance_partition(beta_b$dist_full, beta_b$st, label = "Bacterial")
save_tsv(vp_b$single, "Bact_Variance_single_factor.tsv")
save_tsv(vp_b$multi,  "Bact_Variance_multi_factor.tsv")
save_tsv(vp_b$cohort, "Bact_Variance_cohort_specific.tsv")


# ----------------------------------------------------------------------------
# Step B5  细菌组间距离
# ----------------------------------------------------------------------------

all_distances_bact <- make_all_distances(beta_b$st, beta_b$dist_full)

dist_summary_bact <- all_distances_bact[, .(N = .N,
                                            Mean_Distance   = round(mean(Distance),   4),
                                            SD_Distance     = round(sd(Distance),     4),
                                            Median_Distance = round(median(Distance), 4)),
                                        by = .(Comparison, Mother_Group)]
print(dist_summary_bact)
save_tsv(dist_summary_bact,  "Bact_Inter_group_distances_summary.tsv")
save_tsv(all_distances_bact, "Bact_Inter_group_distances_all.tsv")

p_dist_bact <- plot_dist_box(all_distances_bact, DIST_COMPARISONS,
                             "Bacterial Community Distance: Maternal Groups to Calf Niches")
save_pdf(p_dist_bact, "Fig1_Bact_Inter_Group_Distances.pdf", 8, 4)

dist_tests_bact <- rbindlist(lapply(split(all_distances_bact, all_distances_bact$Comparison), function(d) {
  if (d[, uniqueN(Mother_Group)] < 2) return(NULL)
  kw  <- kruskal.test(Distance ~ Mother_Group, data = d)
  out <- list(data.table(Test = "Kruskal-Wallis",
                         Statistic = round(unname(kw$statistic), 3),
                         P_value = kw$p.value))
  if (kw$p.value < 0.05) {
    pw <- pairwise.wilcox.test(d$Distance, d$Mother_Group, p.adjust.method = "fdr")$p.value
    ij <- which(!is.na(pw), arr.ind = TRUE)
    out[[2]] <- data.table(Test = sprintf("Wilcoxon: %s vs %s",
                                          rownames(pw)[ij[, 1]], colnames(pw)[ij[, 2]]),
                           Statistic = NA_real_, P_value = pw[ij])
  }
  rbindlist(out, fill = TRUE)
}), idcol = "Comparison", fill = TRUE)
print(dist_tests_bact)
save_tsv(dist_tests_bact, "Bact_Inter_group_distance_tests.tsv")


# ############################################################################
#  PART III.  汇总
# ############################################################################

summary_tbl <- rbindlist(list(
  data.table(Dataset = "Virus",     n_OTU = nrow(rpkm),      n_Sample = ncol(rpkm),
             PERMANOVA_R2 = beta_v$r2_permanova, PERMANOVA_P = beta_v$p_permanova,
             Betadisper_P = beta_v$p_disp,
             PCo1_var = beta_v$var_exp[1], PCo2_var = beta_v$var_exp[2]),
  data.table(Dataset = "Bacteria",  n_OTU = nrow(rpkm_bact), n_Sample = ncol(rpkm_bact),
             PERMANOVA_R2 = beta_b$r2_permanova, PERMANOVA_P = beta_b$p_permanova,
             Betadisper_P = beta_b$p_disp,
             PCo1_var = beta_b$var_exp[1], PCo2_var = beta_b$var_exp[2])
), fill = TRUE)

print(summary_tbl)
save_tsv(summary_tbl, "Fig1_summary_virus_vs_bacteria.tsv")

# 拼图（可选）
p_combo <- (beta_v$plot | beta_b$plot) +
  plot_layout(guides = "collect") & theme(legend.position = "right")
save_pdf(p_combo, "Fig1_PCoA_Virus_vs_Bacteria.pdf", 16, 7)

message("\n\u2713 All analyses completed. Outputs in: ", OUT_DIR)



message("Generating polyomavirus plot...")

kraken_data <- fread("../data/fig1/kraken_viruse2_kraken2")
target_spes <- "k__Viruses|k__Shotokuvirae|p__Cossaviricota|c__Papovaviricetes|o__Sepolyvirales|f__Polyomaviridae|s__Bos_taurus_polyomavirus_1"

if (target_spes %in% kraken_data$V1) {
  spes_row <- kraken_data[V1 == target_spes]
  
  poly_data <- data.frame(
    Sample_ID = colnames(kraken_data)[-1],
    Abundance = as.numeric(unlist(spes_row[, -1, with = FALSE])),
    check.names = FALSE
  )
  
  poly_merged <- merge(poly_data, meta[, .(Sample_ID, Group, Area)], by = "Sample_ID")
  poly_merged <- poly_merged[!is.na(poly_merged$Group),]
  
  poly_merged$Group_Area <- paste(poly_merged$Area,poly_merged$Group,sep="_")
  poly_merged$Group_Area <- factor(poly_merged$Group_Area)
  
  poly_colors <- group_color[levels(poly_merged$Group_Area)]
  
  p_poly <- ggplot(poly_merged, aes(x = Group_Area, y = Abundance, fill = Group_Area)) +
    geom_boxplot(outliers = FALSE, linewidth = 0.5) +
    geom_jitter(position = position_jitter(width = 0.15), size = 0.8, color = "#454C49", alpha = 0.7, shape = 20) +
    facet_grid(. ~ Area, scales = "free", space = "free") +
    scale_fill_manual(values = poly_colors) +
    labs(x = NULL, y = "Relative abundance") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
      axis.text.y = element_text(size = 7),
      axis.title.y = element_text(size = 8),
      legend.position = "none",
      strip.text = element_text(size = 8, face = "bold"),
      panel.grid = element_blank()
    )
  
  ggsave("../output/Fig1/Fig1_Polyomavirus.pdf", p_poly, width = 8, height = 5, units = "cm", dpi = 300)
  message("✓ Polyomavirus plot saved")
} else {
  message("⚠ Target species not found in kraken data. Skipping polyomavirus plot.")
}


supp <- do.call(rbind, lapply(unique(poly_merged$Area), function(a) {
  d <- poly_merged[poly_merged$Area == a, ]
  
  do.call(rbind, lapply(c("C0dOR","C0dME"), function(g) {
    do.call(rbind, lapply(c("M0dF","M0dR","M21dF","M21dR","C70dF"), function(m) {
      x <- d$Abundance[d$Group == g]
      y <- d$Abundance[d$Group == m]
      if (!length(x) || !length(y)) return(NULL)
      
      data.frame(
        Area = a, Comparison = paste(g, "vs", m),
        n1 = length(x), n2 = length(y),
        Mean1 = mean(x), Mean2 = mean(y),
        FoldChange = (mean(x)+1e-10)/(mean(y)+1e-10),
        P_value = wilcox.test(x, y, exact=FALSE)$p.value
      )
    }))
  }))
}))

supp$FDR <- p.adjust(supp$P_value, "BH")

write.csv(supp, "../output/Fig1/Supplementary_Polyomavirus.csv", row.names = FALSE)




# ============================================================================
# Step 9  组合图（patchwork）
# ============================================================================

# ----------------------------------------------------------------------------
# 9.1  PCoA 竖排：上 = 病毒，下 = 细菌
# ----------------------------------------------------------------------------
p_pcoa_combo <- (beta_v$plot / beta_b$plot) +
  plot_layout(nrow = 2, guides = "collect") &
  theme(legend.position = "right")

save_pdf(p_pcoa_combo, "Fig1_PCoA_Virus_vs_Bacteria.pdf", 8, 14)

# ----------------------------------------------------------------------------
# 9.2  Shannon 竖排：上 = 病毒，下 = 细菌
# ----------------------------------------------------------------------------
p_shannon_combo <- (p_alpha / p_alpha_bact) +
  plot_layout(nrow = 2) &
  theme(legend.position = "none")

save_pdf(p_shannon_combo, "Fig1_Shannon_Virus_vs_Bacteria.pdf", 7, 9)



p_dist_combo <- (p_dist_viruse / p_dist_bact) +
  plot_layout(nrow = 2) &
  theme(legend.position = "none")
save_pdf(p_dist_combo, "Fig1_dist_Virus_vs_Bacteria.pdf", 7,8)
# ----------------------------------------------------------------------------
# 9.3  四宫格：
#      左上 = Polyomavirus   右上 = CheckV 分布
#      左下 = Step4 丰度     右下 = Step5 病毒科-宿主门
# ----------------------------------------------------------------------------
if (is.null(p_poly)) {
  # 若 polyomavirus 未找到，用占位图避免拼接失败
  p_poly <- ggplot() + annotate("text", 0.5, 0.5,
                                label = "Polyomavirus\nnot detected", size = 4) +
    theme_void()
}

p_panel <- (p_poly   | p_checkv) /
  (p_abund  | p_host) +
  plot_layout(widths = c(1, 1), heights = c(1, 1.2)) +
  plot_annotation(
    tag_levels = list(c("A", "B", "C", "D")),
    theme = theme(plot.tag = element_text(face = "bold", size = 16))
  )

save_pdf(p_panel, "Fig1_Panel_2x2.pdf", 16, 15)
message("\n\u2713 Figure 1 completed successfully")








# ============================================================================
# 快速生成 genus RPTM（= calc_rptm + agg_to_genus 的等价实现）
# ============================================================================
agg_genus <- function(mat, anno_dt, rank) {
  tax <- trimws(as.character(anno_dt[match(rownames(mat), OTU)][[rank]]))
  tax[is.na(tax) | !nzchar(tax)] <- "Unclassified"
  m <- rowsum(mat, group = tax, reorder = TRUE)          # 同属 OTU 按样本求和
  m <- m[order(rowSums(m), decreasing = TRUE), , drop = FALSE]
  data.table(genus = rownames(m), as.data.table(m))      # 首列 genus，兼容下游
}

# 病毒 → 宿主属（对应原脚本 host_annotation2）
anno_NA <- anno[is.na(anno$host_Family)==F&anno$host_Genus!="Unclassified",]
anno2 <- anno_v[anno_v$host_annotation2%in%anno_NA$host_Genus,]
anno2$OTU<- anno2$Contig_name
vir_rptm_genus  <- agg_genus(rpkm,      anno2,      "host_annotation2")

# anno2 <- anno[anno$host_Genus%in%anno_NA$host_Genus,]
# vir_rptm_genus  <- agg_genus(rpkm,      anno2,      "host_Genus")
# 细菌 → 自身属
blast_bact$OTU <- blast_bact$Contig_name
bact_rptm_genus <- agg_genus(rpkm_bact, blast_bact, "annotation_pre")

write.table(vir_rptm_genus,  "../data/fig2_3/vir_rptm_genus.tsv",row.names  = F)
write.table(bact_rptm_genus, "../data/fig2_3/bact_rptm_genus.tsv",row.names = F)

# 校验：属表列和 应等于 OTU 表列和
all.equal(colSums(as.matrix(vir_rptm_genus[, -1])),  colSums(rpkm),      tolerance = 1e-6)
all.equal(colSums(as.matrix(bact_rptm_genus[, -1])), colSums(rpkm_bact), tolerance = 1e-6)

message(sprintf("病毒属: %d x %d | 细菌属: %d x %d",
                nrow(vir_rptm_genus) - 1, ncol(vir_rptm_genus) - 1,
                nrow(bact_rptm_genus) - 1, ncol(bact_rptm_genus) - 1))
