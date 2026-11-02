# ============================================================
# Fig2/Fig2_Temporal_Dynamics.R - Figure 2: Temporal Dynamics (Gansu Cohort)
# Description: LEfSe, heatmaps, Venn diagrams, trend plots, 
#              correlations, pathway differences
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
  library(data.table)
})

source("../R/00_config.R")
source("../R/01_utils.R")
dir.create("../output/Fig2", recursive = TRUE, showWarnings = FALSE)

# ==================== 1. 数据加载 ====================
message("Loading data...")
groups <- c("Gansu_C0dME", "Gansu_C70dF", "Gansu_M0dF")

# 使用 meta.tsv 替代 sample_info.csv
meta <- fread("../data/meta.tsv")
sample_info <- meta[Group_Area %in% groups, 
                    .(Sample_ID, Group_Area)]
rownames(sample_info) <- sample_info$Sample_ID

virus_exp <- fread("../data/fig2_3/vir_rptm_genus.tsv",check.names = F)
virus_exp <- virus_exp[virus_exp$genus!="Unclassified",]
virus_exp_table <- as.data.frame(virus_exp[, -1])
rownames(virus_exp_table) <- virus_exp[[1]]
virus_exp_table <- virus_exp_table[, colnames(virus_exp_table) %in% sample_info$Sample_ID, drop = FALSE]
virus_genus_table <- data.frame(row.names = rownames(virus_exp_table), Genus = rownames(virus_exp_table))

bact_exp <- fread("../data/fig2_3/bact_rptm_genus.tsv",check.names = F)
bact_exp <- bact_exp[bact_exp$genus!="Unclassified",]
bact_exp_table <- as.data.frame(bact_exp[, -1])
rownames(bact_exp_table) <- bact_exp[[1]]
bact_exp_table <- bact_exp_table[, colnames(bact_exp_table) %in% sample_info$Sample_ID, drop = FALSE]
bact_genus_table <- data.frame(row.names = rownames(bact_exp_table), Genus = rownames(bact_exp_table))

spacer <- read.table("../data/fig2_3/spacer_rel_genus.txt",check.names = F)
spacer_df <- as.data.frame(spacer[, -1])
rownames(spacer_df) <- spacer[[1]]
spacer_exp_table <- spacer_df[, colnames(spacer_df) %in% sample_info$Sample_ID, drop = FALSE]

# ==================== 2. LEfSe 分析 ====================
message("Running LEfSe analysis (Gansu cohort)...")

mt_virus <- microtable$new(
  sample_table = data.frame(group = sample_info$Group_Area, row.names = sample_info$Sample_ID),
  otu_table = virus_exp_table,
  tax_table = virus_genus_table
)
set.seed(1234)
t1_virus <- trans_diff$new(dataset = mt_virus, method = "lefse", group = "group", taxa_level = "Genus")
diff_res_virus <- as.data.frame(t1_virus$res_diff)
diff_res_sig_virus <- diff_res_virus[diff_res_virus$LDA > 3 & diff_res_virus$P.adj < 0.01, ]
write.csv(diff_res_sig_virus, "../output/Fig2/Gansu_virus_lefse_genus.csv", row.names = FALSE)

mt_bact <- microtable$new(
  sample_table = data.frame(group = sample_info$Group_Area, row.names = sample_info$Sample_ID),
  otu_table = bact_exp_table,
  tax_table = bact_genus_table
)
set.seed(1234)
t1_bact <- trans_diff$new(dataset = mt_bact, method = "lefse", group = "group", taxa_level = "Genus")
diff_res_bact <- as.data.frame(t1_bact$res_diff)
diff_res_sig_bact <- diff_res_bact[diff_res_bact$LDA > 3 & diff_res_bact$P.adj < 0.01, ]
write.csv(diff_res_sig_bact, "../output/Fig2/Gansu_bact_lefse_genus.csv", row.names = FALSE)

virus_num <- as.data.frame(table(diff_res_sig_virus$Group))
colnames(virus_num) <- c("group", "num")
bact_num <- as.data.frame(table(diff_res_sig_bact$Group))
colnames(bact_num) <- c("group", "num")
total_num <- rbind(virus_num, bact_num)
total_num$type <- c(rep("virus", 3), rep("bacteria", 3))
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
sample_annotation <- sample_annotation[colnames(mat),,drop=F]
virus_annotation <- virus_annotation[rownames(mat),,drop=F]
ann_colors <- list(group = group_color,source=group_color)

ht_virus <- pheatmap(mat, scale="row", show_colnames=F, show_rownames=F,
                     cutree_rows=3, cutree_cols=3,breaks = breaks,
                     annotation_col=sample_annotation, annotation_row=virus_annotation,
                     annotation_colors=ann_colors,
                     color=colorRampPalette(color_palette)(50),
                     cluster_rows=T, cluster_cols=T,
                     clustering_distance_rows="correlation", clustering_distance_cols="correlation",
                     clustering_method="complete",
                     fontsize=10, fontsize_row=8, fontsize_col=10,
                     cellheight=4, cellwidth=4, border_color="white", angle_col="45")

pdf("../output/Fig2/Fig2_Heatmap_virus.pdf", width = 8, height = 8)
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
                    border_color = "white", angle_col = "45")

pdf("../output/Fig2/Fig2_Heatmap_bact.pdf", width = 8, height = 9)
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

# 使用 Gansu_* 格式
birth_virus <- diff_res_sig_virus[diff_res_sig_virus$Group == "Gansu_C0dME", ]$Taxa
birth_bact <- diff_res_sig_bact[diff_res_sig_bact$Group == "Gansu_C0dME", ]$Taxa
plot_venn(birth_virus, birth_bact, "../output/Fig2/Fig2_Venn_birth.pdf")

D70_virus <- diff_res_sig_virus[diff_res_sig_virus$Group == "Gansu_C70dF", ]$Taxa
D70_bact <- diff_res_sig_bact[diff_res_sig_bact$Group == "Gansu_C70dF", ]$Taxa
plot_venn(D70_virus, D70_bact, "../output/Fig2/Fig2_Venn_D70.pdf")

adult_virus <- diff_res_sig_virus[diff_res_sig_virus$Group == "Gansu_M0dF", ]$Taxa
adult_bact <- diff_res_sig_bact[diff_res_sig_bact$Group == "Gansu_M0dF", ]$Taxa
plot_venn(adult_virus, adult_bact, "../output/Fig2/Fig2_Venn_adult.pdf")
message("✓ Venn diagrams saved")

library(VennDiagram); library(gridExtra)
dir.create("../output/Fig2", recursive = TRUE, showWarnings = FALSE)

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
write.table(summary_dt,"../output/Fig2/LEFse_overlap.tsv",quote = F,col.names = F)

# ==================== 5. 趋势图 ====================
message("Generating trend plots...")

# 使用 Gansu_* 作为 group_levels
gansu_levels <- c("Gansu_C0dME", "Gansu_C70dF", "Gansu_M0dF")

p1 <- plot_target("Corynebacterium", virus_exp_table, bact_exp_table, spacer_exp_table,
                  sample_info, "#FF9B80", gansu_levels, gansu_levels)
p2 <- plot_target("Blautia", virus_exp_table, bact_exp_table, spacer_exp_table,
                  sample_info, "#A7C3D9", gansu_levels, gansu_levels)
p3 <- plot_target("Faecalibacterium", virus_exp_table, bact_exp_table, spacer_exp_table,
                  sample_info, "#A7C3D9", gansu_levels, gansu_levels)
p4 <- plot_target("Mediterraneibacter", virus_exp_table, bact_exp_table, spacer_exp_table,
                  sample_info, "#A7C3D9", gansu_levels, gansu_levels)

p_all <- wrap_plots(p2, p3, p4, p1, ncol = 2)
ggsave("../output/Fig2/Fig2_Trend_spacer.pdf", p_all, dpi = 300, width = 10, height = 10, units = "cm")
message("✓ Trend plots saved")

# ==================== 6. 相关性散点图 ====================
message("Generating correlation scatter plots...")

gansu_colors <- c("Gansu_C0dME" = "#FF9B80", "Gansu_C70dF" = "#A7C3D9", "Gansu_M0dF" = "#75CED9")

P1 <- plot_dot_target("Chordicoccus", virus_exp_table, bact_exp_table, sample_info,
                      gansu_levels, gansu_levels, gansu_colors)
P2 <- plot_dot_target("Glutamicibacter", virus_exp_table, bact_exp_table, sample_info,
                      gansu_levels, gansu_levels, gansu_colors)
P3 <- plot_dot_target("Microbacterium", virus_exp_table, bact_exp_table, sample_info,
                      gansu_levels, gansu_levels, gansu_colors)
P4 <- plot_dot_target("Psychrobacter", virus_exp_table, bact_exp_table, sample_info,
                      gansu_levels, gansu_levels, gansu_colors)

P_all <- wrap_plots(P1, P2, P3, P4, ncol = 1)
ggsave("../output/Fig2/Fig2_Correlation_scatter.pdf", P_all, dpi = 300, width = 15, height = 20, units = "cm")
message("✓ Correlation scatter plots saved")

# ==================== 7. 相关性气泡图 ====================
message("Generating correlation bubble plot...")

intersect_merge_list <- strsplit(summary_dt$Shared_Taxa, ";")
intersect_merge_list <- lapply(intersect_merge_list, function(x) gsub(" ", "", x))
names(intersect_merge_list) <- gansu_levels  # 假设行顺序与 baoji_levels 对应

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
  gansu_levels, gansu_levels, taxa_order, enrichment_map
)

print(P_bubble)
ggsave("../output/Fig2/Fig2_Correlation_bubble.pdf", P_bubble, dpi = 300, width = 10, height = 10, units = "cm")
message("✓ Correlation bubble plot saved")


library(Hmisc); library(vegan); library(phyloseq); library(microeco)
library(ggplot2); library(dplyr); library(ggrepel); library(igraph)
library(ggraph); library(RColorBrewer)
library(WGCNA);library(impute);library(preprocessCore)

# ============================================================
# Fig2_Network_Analysis.R - Microbial Interaction Networks
# Description: Co-occurrence networks for Bacteria and Viruses
#              across Gansu cohort groups.
# ============================================================

suppressPackageStartupMessages({
  library(igraph)
  library(microeco)
  library(ggplot2)
  library(ggrepel)
  library(Hmisc) # for rcorr
  library(dplyr)
  library(data.table)
})

source("../R/00_config.R") # 确保加载了 group_color 等全局变量
source("../R/01_utils.R")

# 创建输出目录
net_dir <- "../output/Fig2/Network"
dir.create(net_dir, recursive = TRUE, showWarnings = FALSE)

# ==================== 1. 数据准备 ====================
message("Loading data for Network Analysis...")

# 1.1 加载元数据
meta <- fread("../data/meta.tsv")
groups <- c("Gansu_C0dME", "Gansu_C70dF", "Gansu_M0dF")
sample_info <- meta[Group_Area %in% groups, .(Sample_ID, Group_Area)]
rownames(sample_info) <- sample_info$Sample_ID

# 1.2 加载细菌数据
bact_exp <- fread("../data/fig2_3/bact_rptm_genus.tsv", check.names = F)
bact_exp <- bact_exp[bact_exp$genus != "Unclassified", ]
bact_table <- as.data.frame(bact_exp[, -1])
rownames(bact_table) <- bact_exp[[1]]
# 过滤只保留 sample_info 中存在的样本
bact_table <- bact_table[, colnames(bact_table) %in% sample_info$Sample_ID, drop = FALSE]
bact_taxa <- data.frame(Genus = rownames(bact_table), row.names = rownames(bact_table))

# 1.3 加载病毒数据
virus_exp <- fread("../data/fig2_3/vir_rptm_genus.tsv", check.names = F)
virus_exp <- virus_exp[virus_exp$genus != "Unclassified", ]
virus_table <- as.data.frame(virus_exp[, -1])
rownames(virus_table) <- virus_exp[[1]]
virus_table <- virus_table[, colnames(virus_table) %in% sample_info$Sample_ID, drop = FALSE]
virus_taxa <- data.frame(Genus = rownames(virus_table), row.names = rownames(virus_table))

# ==================== 2. 网络分析核心函数 ====================

# 定义模块颜色映射
global_module_colors <- c(
  "Module_1"="#7FBCB8", "Module_2"="#E6866D", "Module_3"="#A6D854", 
  "Module_4"="#D9A13C", "Module_5"="#B28BB5", "Module_6"="#E3C77B", 
  "Module_7"="#7C9EB2", "Module_8"="#D18F9E", "Module_9"="#A8C9A5", 
  "Module_10"="#E5B375","Module_11"="#7570B3","Module_12"="#99C2A2", 
  "Module_13"="#D9A8B4","Module_14"="#8FB4C3","Module_15"="#C5D68A",
  "Others"="grey70", "Unassigned"="grey50"
)

run_network_analysis <- function(data_table, taxa_table, sample_info, group_name, prefix, output_dir) {
  cat("\n▶ 开始分析:", prefix, "-", group_name, "\n")
  
  # --- 2.1 数据提取 ---
  group_samples <- rownames(sample_info)[sample_info$Group_Area == group_name]
  valid_samples <- intersect(colnames(data_table), group_samples)
  
  if(length(valid_samples) < 3) {
    cat("❌ 跳过:", group_name, "- 有效样本数不足3个 (当前:", length(valid_samples), ")\n")
    return(NULL)
  }
  
  # 提取数据并强制转换为 data.frame (关键修复点)
  genus_df <- as.data.frame(data_table[, valid_samples, drop = FALSE])
  
  # 过滤低丰度 (阈值 0.01%)
  keep_otu <- rowMeans(genus_df) >= 0.0001
  genus_df <- genus_df[keep_otu, , drop = FALSE]
  
  if(nrow(genus_df) < 5) {
    cat("❌ 跳过:", group_name, "- 过滤后OTU过少\n")
    return(NULL)
  }
  cat("  数据维度:", nrow(genus_df), "OTUs x", ncol(genus_df), "Samples\n")
  
  # 准备对应的 tax_table，确保行名匹配
  current_taxa <- taxa_table[rownames(genus_df), , drop = FALSE]
  # 再次确认为 data.frame
  current_taxa <- as.data.frame(current_taxa)
  
  # --- 2.2 相关性计算 (Spearman) ---
  # 将 df 转为 matrix 用于 rcorr
  genus_mat <- as.matrix(genus_df)
  corr_res <- rcorr(t(genus_mat), type = "spearman")
  r_mat <- corr_res$r
  p_mat <- corr_res$P
  
  # 多重检验校正
  p_adj_vec <- p.adjust(as.vector(p_mat), method = "BH")
  p_adj_mat <- matrix(p_adj_vec, nrow = nrow(p_mat))
  dimnames(p_adj_mat) <- dimnames(p_mat)
  
  # 筛选阈值: |r| >= 0.8 & P.adj < 0.01
  adj_mat <- ifelse(abs(r_mat) >= 0.8 & p_adj_mat < 0.01, r_mat, 0)
  diag(adj_mat) <- 0
  
  pos_edges <- sum(adj_mat > 0) / 2
  neg_edges <- sum(adj_mat < 0) / 2
  cat("  边统计: 正相关 =", pos_edges, ", 负相关 =", neg_edges, "\n")
  
  if(pos_edges + neg_edges == 0) {
    cat("❌ 跳过: 无显著相关边\n")
    return(NULL)
  }
  
  # --- 2.3 构建 igraph 对象 ---
  g <- graph_from_adjacency_matrix(adj_mat, weighted = TRUE, mode = "undirected")
  g <- simplify(g) 
  g <- delete_vertices(g, degree(g) == 0) 
  
  if(vcount(g) == 0) return(NULL)
  
  E(g)$correlation <- E(g)$weight
  E(g)$weight <- abs(E(g)$weight)
  
  # 添加物种名称属性
  V(g)$name_orig <- V(g)$name
  V(g)$Genus <- current_taxa[V(g)$name, "Genus"]
  
  # --- 2.4 模块划分 (Fast Greedy) ---
  comm <- cluster_fast_greedy(g, weights = E(g)$weight)
  V(g)$module <- membership(comm)
  V(g)$module_name <- paste0("Module_", V(g)$module)
  
  # 拓扑参数统计
  topo_params <- data.frame(
    Group = group_name, 
    Type = prefix,
    Nodes = vcount(g), 
    Edges = ecount(g),
    Pos_Cor = pos_edges, 
    Neg_Cor = neg_edges,
    Avg_Degree = mean(degree(g)), 
    Density = edge_density(g),
    Modularity = modularity(comm), 
    stringsAsFactors = FALSE
  )
  
  # --- 2.5 Zi-Pi 计算 (使用 microeco) ---
  # 关键修复：确保传入的是 data.frame
  tryCatch({
    dataset_net <- microtable$new(
      otu_table = genus_df, 
      tax_table = current_taxa
    )
    
    t1 <- trans_network$new(dataset = dataset_net, use_WGCNA_pearson_spearman = TRUE,
                            filter_thres = 0.0001, cor_method = "spearman")
    
    # 手动替换为我们已经计算好的 igraph 对象
    t1$res_network <- g
    t1$cal_network_attr()
    t1$cal_module()
    
    net_top <- as.data.frame(t1$get_node_table(node_roles = TRUE)$res_node_table)
    net_top$Genus <- V(g)$Genus[match(net_top$name, V(g)$name_orig)]
    
  }, error = function(e) {
    cat("⚠️ microeco 计算 Zi-Pi 失败，使用基础指标替代:", conditionMessage(e), "\n")
    # 如果 microeco 失败，生成一个基础的 node_table
    net_top <- data.frame(
      name = V(g)$name_orig,
      degree = degree(g),
      betweenness = betweenness(g),
      closeness = closeness(g),
      module = V(g)$module,
      z = NA, p = NA, taxa_roles = "Unknown",
      Genus = V(g)$Genus,
      stringsAsFactors = FALSE
    )
  })
  
  # --- 2.6 可视化准备 ---
  mod_counts <- sort(table(V(g)$module_name), decreasing = TRUE)
  top_mods <- names(mod_counts)[1:min(10, length(mod_counts))]
  
  V(g)$mod_simple <- ifelse(V(g)$module_name %in% top_mods, V(g)$module_name, "Others")
  V(g)$color <- global_module_colors[V(g)$mod_simple]
  V(g)$color[is.na(V(g)$color)] <- "grey50"
  
  E(g)$edge_color <- ifelse(E(g)$correlation > 0, "#FF6B6B", "#4ECDC4")
  
  deg <- degree(g)
  if(max(deg) > min(deg)) {
    V(g)$size <- 3 * sqrt(1 + 3 * (deg - min(deg)) / max(1, max(deg) - min(deg)))
  } else {
    V(g)$size <- 3
  }
  
  # --- 2.7 保存结果 ---
  
  # A. 保存拓扑参数
  write.csv(topo_params, file.path(output_dir, paste0(prefix, "_topology_", group_name, ".csv")), row.names = FALSE)
  
  # B. 保存节点属性表 (Zi-Pi)
  write.csv(net_top, file.path(output_dir, paste0(prefix, "_node_ZiPi_", group_name, ".csv")), row.names = FALSE)
  
  # C. 绘制网络图 (PDF)
  pdf_file <- file.path(output_dir, paste0(prefix, "_network_", group_name, ".pdf"))
  pdf(pdf_file, width = 12, height = 12)
  plot(g, 
       layout = layout_with_fr(g), 
       vertex.label = NA, 
       vertex.frame.color = "white",
       vertex.color = V(g)$color,
       vertex.size = V(g)$size,
       edge.width = 0.5, 
       edge.curved = 0.1, 
       edge.color = E(g)$edge_color,
       main = paste(prefix, group_name, sep = " - "))
  legend("topright", legend = c("Positive Cor", "Negative Cor"), 
         col = c("#FF6B6B","#4ECDC4"), lwd = 2, bty = "n", cex = 0.8)
  dev.off()
  
  # D. 绘制 Zi-Pi 散点图 (PDF) - ✅ 修复颜色问题
  zipi_file <- file.path(output_dir, paste0(prefix, "_ZiPi_", group_name, ".pdf"))
  
  if("z" %in% colnames(net_top) && "p" %in% colnames(net_top)) {
    # 1. 创建简化的模块标签 (M1, M2...)
    net_top$mod_simple_plot <- factor(
      ifelse(net_top$module %in% gsub("Module_","M",top_mods), gsub("Module_","M",net_top$module), "Others"),
      levels = c(gsub("Module_","M",top_mods), "Others")
    )
    
    # 2. ✅ 关键修复：创建与简化标签匹配的颜色向量
    zipi_colors <- global_module_colors
    names(zipi_colors) <- gsub("Module_", "M", names(global_module_colors))
    # 确保 "Others" 也有颜色
    if(!"Others" %in% names(zipi_colors)) zipi_colors["Others"] <- "grey70"
    
    p_zipi <- ggplot(net_top, aes(x = p, y = z, color = mod_simple_plot)) +
      geom_vline(xintercept = 0.62, linetype = "dashed", color = "grey50") +
      geom_hline(yintercept = 2.5, linetype = "dashed", color = "grey50") +
      geom_point(size = 1.5, alpha = 0.8) +
      # 使用修复后的 zipi_colors
      scale_color_manual(values = zipi_colors, drop = FALSE, name = "Modules") +
      geom_text_repel(data = subset(net_top, !is.na(z) & !is.na(p) & (z > 2.5 | p > 0.62)),
                      aes(label = Genus), size = 2.5, max.overlaps = 20, show.legend = FALSE) +
      labs(x = "Among-module connectivity (Pi)", y = "Within-module connectivity (Zi)",
           title = paste(prefix, group_name)) +
      theme_classic() + 
      theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 10)) +
      xlim(-0.1, 0.8) + ylim(-2, 3)
    
    ggsave(zipi_file, p_zipi, width =9, height = 6, units = "cm")
  } else {
    cat("⚠️ 跳过 Zi-Pi 绘图，因为缺少 z/p 值\n")
  }
  
  cat("✅", prefix, group_name, "分析完成!\n")
  
  return(list(graph = g, node_table = net_top, topology = topo_params))
}

# ==================== 3. 批量执行 ====================

all_results <- list()
group_list <- unique(sample_info$Group_Area)

# 3.1 运行细菌网络
cat("\n==================== 细菌网络分析 ====================\n")
for(gname in group_list) {
  res <- tryCatch({
    run_network_analysis(bact_table, bact_taxa, sample_info, gname, "Bact", net_dir)
  }, error = function(e) {
    cat("❌ Bact", gname, "报错:", conditionMessage(e), "\n")
    NULL
  })
  if(!is.null(res)) all_results[[paste0("Bact_", gname)]] <- res
}

# 3.2 运行病毒网络
cat("\n==================== 病毒网络分析 ====================\n")
for(gname in group_list) {
  res <- tryCatch({
    run_network_analysis(virus_table, virus_taxa, sample_info, gname, "Virus", net_dir)
  }, error = function(e) {
    cat("❌ Virus", gname, "报错:", conditionMessage(e), "\n")
    NULL
  })
  if(!is.null(res)) all_results[[paste0("Virus_", gname)]] <- res
}

cat("\n🎉 所有网络分析已完成！结果保存在:", net_dir, "\n")

