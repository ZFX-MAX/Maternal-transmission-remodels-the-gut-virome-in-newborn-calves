# ============================================================
# Fig3/Fig3_Niche_Differences.R - Figure 3: Niche Differences (Baoji Cohort)
# Description: LEfSe, heatmaps, Venn diagrams, Procrustes,
#              diversity correlation, trend plots, correlations
# ============================================================
suppressPackageStartupMessages({
  library(pheatmap)
  library(gtable)
  library(microeco)
  library(ggplot2)
  library(ggview)
  library(VennDiagram)
  library(dplyr)
  library(tidyr)
  library(grid)
  library(patchwork)
  library(vegan)
  library(psych)
  library(data.table)
})

source("../R/00_config.R")
source("../R/01_utils.R")
dir.create("../output/Fig3", recursive = TRUE, showWarnings = FALSE)

# ==================== 1. 数据加载 ====================
message("Loading data...")

# 使用 meta.tsv 替代 sample_info.csv
meta <- fread("../data/meta.tsv")
sample_info <- meta[Group_Area %in% c("Baoji_C0dOR", "Baoji_C0dME", "Baoji_M0dR", "Baoji_M0dF"), 
                    .(Sample_ID, Group_Area)]
rownames(sample_info) <- sample_info$Sample_ID
group_color <- group_color[names(group_color)%in%sample_info$Group_Area]

virus_exp <- fread("../data/fig2_3/vir_rptm_genus.tsv")
virus_exp_table <- as.data.frame(virus_exp[, -1])
rownames(virus_exp_table) <- virus_exp[[1]]
virus_exp_table <- virus_exp_table[, colnames(virus_exp_table) %in% sample_info$Sample_ID, drop = FALSE]
virus_genus_table <- data.frame(row.names = rownames(virus_exp_table), Genus = rownames(virus_exp_table))

bact_exp <- fread("../data/fig2_3/bact_rptm_genus.tsv")
bact_exp_table <- as.data.frame(bact_exp[, -1])
rownames(bact_exp_table) <- bact_exp[[1]]
bact_exp_table <- bact_exp_table[, colnames(bact_exp_table) %in% sample_info$Sample_ID, drop = FALSE]
bact_genus_table <- data.frame(row.names = rownames(bact_exp_table), Genus = rownames(bact_exp_table))

spacer <- read.table("../data/fig2_3/spacer_rel_genus.txt",check.names = F)
spacer_df <- as.data.frame(spacer[, -1])
rownames(spacer_df) <- spacer[[1]]
spacer_exp_table <- spacer_df[, colnames(spacer_df) %in% sample_info$Sample_ID, drop = FALSE]

# ==================== 2. LEfSe 分析 ====================
message("Running LEfSe analysis (Baoji cohort)...")

mt_virus <- microtable$new(
  sample_table = data.frame(group = sample_info$Group_Area, row.names = sample_info$Sample_ID),
  otu_table = virus_exp_table,
  tax_table = virus_genus_table
)
set.seed(1234)
t1_virus <- trans_diff$new(dataset = mt_virus, method = "lefse", group = "group", taxa_level = "Genus")
diff_res_virus <- as.data.frame(t1_virus$res_diff)
diff_res_sig_virus <- diff_res_virus[diff_res_virus$LDA > 3 & diff_res_virus$P.adj < 0.01, ]
write.csv(diff_res_sig_virus, "../output/Fig3/Baoji_virus_lefse_genus.csv", row.names = FALSE)

mt_bact <- microtable$new(
  sample_table = data.frame(group = sample_info$Group_Area, row.names = sample_info$Sample_ID),
  otu_table = bact_exp_table,
  tax_table = bact_genus_table
)
set.seed(1234)
t1_bact <- trans_diff$new(dataset = mt_bact, method = "lefse", group = "group", taxa_level = "Genus")
diff_res_bact <- as.data.frame(t1_bact$res_diff)
diff_res_sig_bact <- diff_res_bact[diff_res_bact$LDA > 3 & diff_res_bact$P.adj < 0.01, ]
write.csv(diff_res_sig_bact, "../output/Fig3/Baoji_bact_lefse_genus.csv", row.names = FALSE)

virus_num <- as.data.frame(table(diff_res_sig_virus$Group))
colnames(virus_num) <- c("group", "num")
bact_num <- as.data.frame(table(diff_res_sig_bact$Group))
colnames(bact_num) <- c("group", "num")
total_num <- rbind(virus_num, bact_num)
total_num$type <- c(rep("virus", 4), rep("bacteria", 4))
message("✓ LEfSe analysis completed")

# ==================== 3. 热图 ====================
message("Generating heatmaps...")

highlight_genes <- c("Acinetobacter", "Corynebacterium", "Glutamicibacter",
                     "Microbacterium", "Psychrobacter", "Aerococcus", "Bifidobacterium")
color_palette <- c("#2166ac", "#f4f5f4", "#b2182b")
breaks <- seq(-2, 2, length.out = 50)

# 病毒热图
virus_enriched <- virus_exp_table[rownames(virus_exp_table) %in% diff_res_sig_virus$Taxa, ]

virus_annotation <- data.frame(genus = diff_res_sig_virus$Taxa, source = diff_res_sig_virus$Group)
rownames(virus_annotation) <- virus_annotation$genus
virus_annotation$genus <- NULL

sample_annotation <- data.frame(group = sample_info$Group_Area, row.names = sample_info$Sample_ID)
mat <- as.matrix(virus_enriched)
sample_annotation <- sample_annotation[colnames(mat), , drop = FALSE]
virus_annotation <- virus_annotation[rownames(mat), , drop = FALSE]

ann_colors <- list(group = group_color, source = group_color)

ht_virus <- pheatmap(mat, scale = "row", show_colnames = FALSE, show_rownames = FALSE,
                     cutree_rows = 2, cutree_cols = 2,breaks = breaks,
                     annotation_col = sample_annotation, annotation_row = virus_annotation,
                     annotation_colors = ann_colors,
                     color = colorRampPalette(color_palette)(50),
                     cluster_rows = TRUE, cluster_cols = TRUE,
                     clustering_distance_rows = "correlation", clustering_distance_cols = "correlation",
                     clustering_method = "complete",
                     fontsize = 10, fontsize_row = 8, fontsize_col = 10,
                     cellheight = 4, cellwidth = 4,
                     border_color = "white", angle_col = "45",annotation_names_row = F)

pdf("../output/Fig3/Fig3_Heatmap_virus.pdf", width = 8, height = 8)
ht_virus
dev.off()


bact_enriched <- bact_exp_table[rownames(bact_exp_table) %in% diff_res_sig_bact$Taxa, ]
bact_annotation <- data.frame(genus = diff_res_sig_bact$Taxa, source = diff_res_sig_bact$Group)
rownames(bact_annotation) <- bact_annotation$genus
bact_annotation$genus <- NULL

sample_annotation <- data.frame(group = sample_info$Group_Area, row.names = sample_info$Sample_ID)
mat <- as.matrix(bact_enriched)
sample_annotation <- sample_annotation[colnames(mat), , drop = FALSE]
bact_annotation <- bact_annotation[rownames(mat), , drop = FALSE]

ann_colors <- list(group = group_color, source = group_color)

ht_bact <- pheatmap(mat, scale = "row", show_colnames = FALSE, show_rownames = FALSE,
                    cutree_rows = 4, cutree_cols = 4,breaks = breaks,
                    annotation_col = sample_annotation, annotation_row = bact_annotation,
                    annotation_colors = ann_colors,
                    color = colorRampPalette(color_palette)(50),
                    cluster_rows = TRUE, cluster_cols = TRUE,
                    clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
                    clustering_method = "ward.D2",
                    fontsize = 10, fontsize_row = 8, fontsize_col = 10,
                    cellheight = 4, cellwidth = 4,
                    border_color = "white", angle_col ="45")

pdf("../output/Fig3/Fig3_Heatmap_bact.pdf", width = 8.5, height = 9)
ht_bact
dev.off()
message("✓ Heatmaps saved")

# ==================== 4. 韦恩图 ====================
message("Generating Venn diagrams...")

plot_venn <- function(virus_taxa, bact_taxa, filename) {
  p <- venn.diagram(
    x = list(Virus = virus_taxa, Bacteria = bact_taxa),
    filename = NULL, scaled = FALSE,
    fill = c("#6e48fb", "#ffa500"), cat.col = c("#6e48fb", "#ffa500"),
    alpha = c(0.5, 0.5), col = "black", lwd = 1.3, lty = 1,
    cex = 1.5, cat.cex = 1.3, cat.default.pos = "outer", cat.dist = c(0.15, 0.15), margin = 0.25
  )
  pdf(filename, width = 4, height = 4)
  grid.newpage()
  grid.draw(p)
  dev.off()
}

# 使用 Baoji_* 格式
BCO_virus <- diff_res_sig_virus[diff_res_sig_virus$Group == "Baoji_C0dOR", ]$Taxa
BCO_bact <- diff_res_sig_bact[diff_res_sig_bact$Group == "Baoji_C0dOR", ]$Taxa
plot_venn(BCO_virus, BCO_bact, "../output/Fig3/Fig3_Venn_BCO.pdf")

BCM_virus <- diff_res_sig_virus[diff_res_sig_virus$Group == "Baoji_C0dME", ]$Taxa
BCM_bact <- diff_res_sig_bact[diff_res_sig_bact$Group == "Baoji_C0dME", ]$Taxa
plot_venn(BCM_virus, BCM_bact, "../output/Fig3/Fig3_Venn_BCM.pdf")

BMF0_virus <- diff_res_sig_virus[diff_res_sig_virus$Group == "Baoji_M0dF", ]$Taxa
BMF0_bact <- diff_res_sig_bact[diff_res_sig_bact$Group == "Baoji_M0dF", ]$Taxa
plot_venn(BMF0_virus, BMF0_bact, "../output/Fig3/Fig3_Venn_BMF0.pdf")

BMR0_virus <- diff_res_sig_virus[diff_res_sig_virus$Group == "Baoji_M0dR", ]$Taxa
BMR0_bact <- diff_res_sig_bact[diff_res_sig_bact$Group == "Baoji_M0dR", ]$Taxa
plot_venn(BMR0_virus, BMR0_bact, "../output/Fig3/Fig3_Venn_BMR0.pdf")
message("✓ Venn diagrams saved")


groups <- c("Baoji_C0dOR", "Baoji_C0dME", "Baoji_M0dF", "Baoji_M0dR")

tx <- function(df, g) unique(df$Taxa[df$Group == g])

summary_dt <- rbindlist(lapply(groups, function(g) {
  v <- tx(diff_res_sig_virus, g)
  b <- tx(diff_res_sig_bact,  g)
  data.table(
    Group       = g,
    Virus_N     = length(v),
    Bact_N      = length(b),
    Shared_N    = length(intersect(v, b)),
    Shared_Taxa = paste(intersect(v, b), collapse = ";")
  )
}))

print(summary_dt)
write.table(summary_dt,"../output/Fig3/LEFse_overlap.tsv",quote = F,col.names = F)
message("Generating Procrustes analysis...")

virus <- fread("../data/fig2_3/vir_rptm_genus.tsv")
virus_df <- as.data.frame(virus[, -1])
rownames(virus_df) <- virus[[1]]

bact <- fread("../data/fig2_3/bact_rptm_genus.tsv")
bact_df <- as.data.frame(bact[, -1])
rownames(bact_df) <- bact[[1]]

# 使用 Group_Area 进行筛选
pro_analesis <- function(group_areas, plot_title) {
  sample_select <- meta[Group_Area %in% group_areas, ]
  virus_filter <- virus_df[, colnames(virus_df) %in% sample_select$Sample_ID]
  bact_filter <- bact_df[, colnames(bact_df) %in% sample_select$Sample_ID]
  
  common_sample <- intersect(colnames(virus_filter), colnames(bact_filter))
  virus_filter <- virus_filter[, common_sample]
  bact_filter <- bact_filter[, common_sample]
  
  virus_filter <- virus_filter[rowSums(virus_filter) > 0, ]
  bact_filter <- bact_filter[rowSums(bact_filter) > 0, ]
  
  virus_filter <- as.data.frame(t(virus_filter))
  bact_filter <- as.data.frame(t(bact_filter))
  
  p <- procrustes_analysis(virus_filter, bact_filter, plot_title)
  return(p)
}

p1 <- pro_analesis(c("Baoji_C0dOR"), "Baoji Oral - virus and bacteria")
p2 <- pro_analesis(c("Baoji_C0dME"), "Baoji Meconium - virus and bacteria")
p3 <- pro_analesis(c("Baoji_M0dF"), "Baoji Feces - virus and bacteria")
p4 <- pro_analesis(c("Baoji_M0dR"), "Baoji Rumen - virus and bacteria")

p_all <- wrap_plots(p1, p2, ncol = 2)
ggsave("../output/Fig3/Fig3_Procrustes.pdf", p_all, dpi = 300, width = 10, height = 5, units = "cm")
message("✓ Procrustes analysis saved")

# ==================== 6. 多样性相关性 ====================
message("Generating diversity correlation plots...")

Shannon <- rbind( data.frame( read.table("../output/Fig1/Bact_alpha_diversity.tsv",header = T),type="Bacteria" ),
                  data.frame( read.table("../output/Fig1/Viral_alpha_diversity.tsv",header = T),type="Viruse" )  )
# colnames(Shannon)[colnames(Shannon)=="shannon"] <- "diversity"
# 使用 Group_Area 格式
diversity_correlation <- function(shannon_data, group_area, color) {
  # 从 Group_Area 中提取 Area 和 Group
  parts <- strsplit(as.character(group_area), "_")[[1]]
  area <- parts[1]
  group <- parts[2]
  
  # 确保数据类型一致
  shannon_data$Area <- as.character(shannon_data$Area)
  shannon_data$Group <- as.character(shannon_data$Group)
  shannon_data$type <- as.character(shannon_data$type)
  shannon_data$Sample_ID <- as.character(shannon_data$Sample_ID)
  
  # 筛选当前组的所有数据
  group_data <- shannon_data[shannon_data$Group_Area==group_area, ]
  
  if (nrow(group_data) == 0) {
    warning(paste("No data found for", group_area))
    return(NULL)
  }
  
  # 找出同时有 Bacteria 和 Viruse 数据的 Sample_ID
  bac_samples <- unique(group_data$Sample_ID[group_data$type == "Bacteria"])
  vir_samples <- unique(group_data$Sample_ID[group_data$type == "Viruse"])
  common_samples <- intersect(bac_samples, vir_samples)
  
  if (length(common_samples) < 3) {
    warning(paste("Insufficient matched samples for", group_area, 
                  "(", length(common_samples), "samples)"))
    return(NULL)
  }
  
  # 提取共同样本的细菌和病毒多样性
  bac_data <- group_data[group_data$Sample_ID %in% common_samples & 
                           group_data$type == "Bacteria", 
                         c("Sample_ID", "shannon")]
  names(bac_data)[2] <- "Bacteria"
  
  vir_data <- group_data[group_data$Sample_ID %in% common_samples & 
                           group_data$type == "Viruse", 
                         c("Sample_ID", "shannon")]
  names(vir_data)[2] <- "Viruse"
  
  # 按 Sample_ID 排序以确保一一对应
  bac_data <- bac_data[order(bac_data$Sample_ID), ]
  vir_data <- vir_data[order(vir_data$Sample_ID), ]
  
  # 直接合并（此时 Sample_ID 应该完全对应）
  merge_data <- merge(bac_data, vir_data, by = "Sample_ID")
  
  if (nrow(merge_data) < 3) {
    warning(paste("Final merged data insufficient for", group_area))
    return(NULL)
  }
  
  # 移除 Sample_ID 列用于绘图
  plot_data <- merge_data[, c("Bacteria", "Viruse")]
  
  # 计算相关性
  cor.result <- corr.test(plot_data, method = "spearman", adjust = "none")
  r <- cor.result$r[1, 2]
  p <- cor.result$p[1, 2]
  
  stat_text <- paste0("r = ", round(r, 4), ", p = ", round(p, 4))
  title_text <- paste(group_area, "Alpha Diversity Correlation")
  
  # 绘图
  p_plot <- ggplot(merge_data, aes(x = Bacteria, y = Viruse)) +
    geom_point(shape = 16, size = 1.5, color = color, alpha = 0.7) +
    stat_smooth(method = lm, formula = y ~ x, se = TRUE, 
                color = "black", size = 0.5, alpha = 0.2) +
    theme_bw() +
    labs(title = title_text, 
         x = "Bacterial Shannon Diversity", 
         y = "Viral Shannon Diversity") +
    theme(
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9),
      axis.line = element_line(color = "black"),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5)
    ) +
    annotate("text", 
             x = min(merge_data$Bacteria, na.rm = TRUE),
             y = max(merge_data$Viruse, na.rm = TRUE),
             label = stat_text, 
             hjust = 0, vjust = 1, 
             size = 3, color = "black", fontface = "italic")
  
  cat("✓", group_area, "-", nrow(merge_data), "matched samples\n")
  
  return(p_plot)
}

p1 <- diversity_correlation(shannon_data = Shannon, group_area = "Baoji_C0dOR",color = group_color[["Baoji_C0dOR"]])
p2 <- diversity_correlation(Shannon, "Baoji_C0dME", group_color[["Baoji_C0dME"]])
p3 <- diversity_correlation(Shannon, "Baoji_M0dR", group_color[["Baoji_M0dR"]])
p4 <- diversity_correlation(Shannon, "Baoji_M0dF", group_color[["Baoji_M0dF"]] )

p_all <- wrap_plots(p1, p2, p3, p4, ncol = 2)
ggsave("../output/Fig3/Fig3_Diversity_Correlation.pdf", p_all, dpi = 300, width = 12, height = 12, units = "cm")
message("✓ Diversity correlation plots saved")

# ==================== 7. 趋势图 ====================
message("Generating trend plots...")

# 使用 Baoji_* 作为 group_levels
baoji_levels <- c("Baoji_C0dOR", "Baoji_C0dME", "Baoji_M0dR", "Baoji_M0dF")

p1 <- plot_target("Staphylococcus", virus_exp_table, bact_exp_table, spacer_exp_table,
                  sample_info, "#FF9B80", baoji_levels, baoji_levels)
p2 <- plot_target("Psychrobacter", virus_exp_table, bact_exp_table, spacer_exp_table,
                  sample_info, "#FF9B80", baoji_levels, baoji_levels)
p3 <- plot_target("Jeotgalicoccus", virus_exp_table, bact_exp_table, spacer_exp_table,
                  sample_info, "#FF9B80", baoji_levels, baoji_levels)
p4 <- plot_target("Glutamicibacter", virus_exp_table, bact_exp_table, spacer_exp_table,
                  sample_info, "#FF9B80", baoji_levels,baoji_levels)
p5 <- plot_target("Corynebacterium", virus_exp_table, bact_exp_table, spacer_exp_table,
                  sample_info, "#FF9B80", baoji_levels, baoji_levels)
p6 <- plot_target("Brachybacterium", virus_exp_table, bact_exp_table, spacer_exp_table,
                  sample_info, "#FF9B80", baoji_levels, baoji_levels)
p7 <- plot_target("Aerococcus", virus_exp_table, bact_exp_table, spacer_exp_table,
                  sample_info, "#FF9B80", baoji_levels, baoji_levels)

p_all <- wrap_plots(p2, p4, p5, p7, ncol = 2)
ggsave("../output/Fig3/Fig3_Trend_spacer_oral.pdf", p_all, dpi = 300, width = 15, height = 15, units = "cm")
message("✓ Trend plots saved")

# ==================== 8. 相关性散点图 ====================
message("Generating correlation scatter plots...")

P1 <- plot_dot_target("Alistipes", virus_exp_table, bact_exp_table, sample_info,
                      baoji_levels, baoji_levels, group_color)
P2 <- plot_dot_target("Glutamicibacter", virus_exp_table, bact_exp_table, sample_info,
                      baoji_levels, baoji_levels, group_color)
P3 <- plot_dot_target("Staphylococcus", virus_exp_table, bact_exp_table, sample_info,
                      baoji_levels, baoji_levels, group_color)
P4 <- plot_dot_target("Psychrobacter", virus_exp_table, bact_exp_table, sample_info,
                      baoji_levels, baoji_levels, group_color)

P_all <- wrap_plots(P1, P2, P3, P4, ncol = 1)
ggsave("../output/Fig3/Fig3_Correlation_scatter.pdf", P_all, dpi = 300, width = 15, height = 20, units = "cm")
message("✓ Correlation scatter plots saved")

# ==================== 9. 相关性气泡图 ====================
message("Generating correlation bubble plot...")
intersect_merge_list <- strsplit(summary_dt$Shared_Taxa, ";")
intersect_merge_list <- lapply(intersect_merge_list, function(x) gsub(" ", "", x))
baoji_levels <- c("Baoji_C0dOR", "Baoji_C0dME", "Baoji_M0dF", "Baoji_M0dR")
names(intersect_merge_list) <- baoji_levels  # 假设行顺序与 baoji_levels 对应

intersect_merge <- unlist(intersect_merge_list)
taxa_order <- intersect_merge
enrichment_map <- data.frame(
  taxa = unlist(intersect_merge_list),
  Enriched_Group = rep(names(intersect_merge_list), 
                       sapply(intersect_merge_list, length)),
  stringsAsFactors = FALSE
)
P_bubble <- plot_correlation_bubble(
  intersect_merge, virus_exp_table, bact_exp_table, sample_info,
  baoji_levels, baoji_levels, taxa_order, enrichment_map
)

print(P_bubble)

ggsave("../output/Fig3/Fig3_Correlation_bubble.pdf", P_bubble, dpi = 300, width = 10, height = 10, units = "cm")
message("✓ Correlation bubble plot saved")

message("✓ Figure 3 completed")


