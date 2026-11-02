# ============================================================
# R/01_utils.R - 公共函数库（更新版）
# Description: 提取 Figure 1-3 中的重复函数
# ============================================================
suppressPackageStartupMessages({
  library(microeco)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(patchwork)
  library(vegan)
  library(psych)
  library(GGally)
})

# ==================== 1. LEfSe 分析函数 ====================
run_lefse_genus <- function(otu_file, sample_ids, group_labels, organism_type, output_dir) {
  message(sprintf("Running LEfSe for %s (genus level)...", organism_type))
  
  otu <- read.csv(otu_file, sep = "\t", check.names = FALSE, row.names = 1)
  otu_table <- otu[, colnames(otu) %in% sample_ids, drop = FALSE]
  
  genus_table <- data.frame(Genus = rownames(otu_table), row.names = rownames(otu_table))
  
  sample_table <- data.frame(
    group = group_labels,
    row.names = sample_ids
  )
  
  mt <- microtable$new(
    sample_table = sample_table,
    otu_table = otu_table,
    tax_table = genus_table
  )
  
  set.seed(1234)
  lefse_res <- trans_diff$new(
    dataset = mt,
    method = "lefse",
    group = "group",
    taxa_level = "Genus"
  )
  
  diff_res <- as.data.frame(lefse_res$res_diff)
  diff_res_sig <- diff_res[diff_res$LDA > 3 & diff_res$P.adj < 0.01, ]
  
  output_file <- file.path(output_dir, sprintf("%s_lefse_genus.csv", organism_type))
  write.csv(diff_res_sig, output_file, row.names = FALSE)
  
  message(sprintf("  Found %d significant taxa", nrow(diff_res_sig)))
  return(list(res = diff_res_sig, otu_table = otu_table, sample_table = sample_table))
}

# ==================== 2. 趋势图函数 ====================
plot_target <- function(target, virus_exp_table, bact_exp_table, spacer_exp_table,
                        sample_info, col, group_levels, stage_labels) {
  target_virus <- virus_exp_table[rownames(virus_exp_table) %in% target, , drop = FALSE]
  target_bact <- bact_exp_table[rownames(bact_exp_table) %in% target, , drop = FALSE]
  target_spacer <- spacer_exp_table[rownames(spacer_exp_table) %in% target, , drop = FALSE]
  
  virus_long <- target_virus %>%
    as.data.frame() %>%
    rownames_to_column("Taxa") %>%
    pivot_longer(cols = -Taxa, names_to = "Sample_ID", values_to = "abundance") %>%
    mutate(abundance = as.numeric(abundance), type = "Virus")
  
  bact_long <- target_bact %>%
    as.data.frame() %>%
    rownames_to_column("Taxa") %>%
    pivot_longer(cols = -Taxa, names_to = "Sample_ID", values_to = "abundance") %>%
    mutate(abundance = as.numeric(abundance), type = "Bacterium")
  
  spacer_long <- target_spacer %>%
    as.data.frame() %>%
    rownames_to_column("Taxa") %>%
    pivot_longer(cols = -Taxa, names_to = "Sample_ID", values_to = "abundance") %>%
    mutate(abundance = as.numeric(abundance), type = "Spacer")
  
  plot_data <- bind_rows(virus_long, bact_long, spacer_long) %>%
    left_join(sample_info %>% select(Sample_ID, Group_Area), by = "Sample_ID") %>%
    filter(!is.na(Group_Area)) %>%
    mutate(
      Group_Area = factor(Group_Area, levels = group_levels, labels = stage_labels),
      type = factor(type, levels = c("Virus", "Bacterium", "Spacer"))
    )
  
  plot_summary <- plot_data %>%
    group_by(Group_Area, type) %>%
    summarise(
      mean_abundance = mean(abundance, na.rm = TRUE),
      sd_abundance = sd(abundance, na.rm = TRUE),
      n = sum(!is.na(abundance)),
      se_abundance = ifelse(n > 1, sd_abundance / sqrt(n), 0),
      .groups = "drop"
    ) %>%
    filter(type == "Spacer")
  
  p <- ggplot(plot_summary, aes(x = Group_Area, y = mean_abundance, group = type,
                                color = type, shape = type)) +
    geom_line(linewidth = 0.5, na.rm = TRUE) +
    geom_point(size = 1, na.rm = TRUE) +
    geom_errorbar(aes(ymin = pmax(mean_abundance - se_abundance, 0),
                      ymax = mean_abundance + se_abundance),
                  width = 0.12, linewidth = 0.6, na.rm = TRUE) +
    scale_color_manual(values = c("Spacer" = col)) +
    labs(title = target, x = NULL, y = "Relative abundance") +
    facet_wrap(~type, scales = "free_y", nrow = 1) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 6, face = "italic", hjust = 0.5),
      axis.text = element_text(size = 6, angle = 45, hjust = 1, vjust = 1),
      axis.title = element_text(size = 7),
      strip.text = element_text(size = 5, face = "bold"),
      legend.position = "none"
    )
  
  return(p)
}

# ==================== 3. 相关性散点图函数 ====================
plot_dot_target <- function(target, virus_exp_table, bact_exp_table, sample_info,
                            group_levels, stage_labels, colors) {
  target_virus <- virus_exp_table[rownames(virus_exp_table) %in% target, , drop = FALSE]
  target_bact <- bact_exp_table[rownames(bact_exp_table) %in% target, , drop = FALSE]
  
  virus_long <- target_virus %>%
    as.data.frame() %>%
    rownames_to_column("Taxa") %>%
    pivot_longer(cols = -Taxa, names_to = "Sample_ID", values_to = "virus") %>%
    mutate(type = "Virus")
  
  bact_long <- target_bact %>%
    as.data.frame() %>%
    rownames_to_column("Taxa") %>%
    pivot_longer(cols = -Taxa, names_to = "Sample_ID", values_to = "bact") %>%
    mutate(type = "Bacterium")
  
  final <- merge(virus_long[, c("Sample_ID", "virus")],
                 bact_long[, c("Sample_ID", "bact")], by = "Sample_ID")
  final_group <- merge(final, sample_info %>% select(Sample_ID, Group_Area), by = "Sample_ID")
  
  correlation_result <- data.frame()
  for (i in seq_along(group_levels)) {
    stage_data <- final_group[final_group$Group_Area == group_levels[i], ]
    if (nrow(stage_data) < 3) next
    
    cor_test <- cor.test(stage_data$virus, stage_data$bact, method = "spearman")
    correlation_result <- rbind(correlation_result, data.frame(
      Stage = stage_labels[i],
      rho = cor_test$estimate,
      p_value = cor_test$p.value
    ))
  }
  
  correlation_result <- correlation_result %>%
    mutate(
      Stage = factor(Stage, levels = stage_labels),
      label = paste0("r = ", sprintf("%.2f", rho), "\nP = ", p_value)
    )
  
  final_group <- final_group %>%
    mutate(Stage = factor(Group_Area, levels = group_levels, labels = stage_labels))
  
  p <- ggplot(final_group, aes(x = virus, y = bact, fill = Stage)) +
    geom_point(shape = 21, size = 2, alpha = 0.8, color = "black", stroke = 0.3) +
    geom_smooth(aes(color = Stage), method = "lm", se = TRUE, linewidth = 0.5, show.legend = FALSE) +
    geom_text(data = correlation_result, aes(x = -Inf, y = Inf, label = label),
              inherit.aes = FALSE, hjust = -0.05, vjust = 1.1, size = 3) +
    facet_wrap(~Stage, nrow = 1, scales = "free") +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    labs(x = "Virus abundance", y = "Bacterial abundance") +
    theme_classic() +
    theme(
      axis.text = element_text(size = 7),
      axis.title = element_text(size = 8),
      strip.text = element_text(size = 7, face = "bold"),
      strip.background = element_blank(),
      legend.position = "none"
    )
  
  return(p)
}

# ==================== 4. 相关性气泡图函数 ====================
plot_correlation_bubble <- function(intersect_merge, virus_exp_table,
                                    bact_exp_table, sample_info,
                                    group_levels, stage_labels, taxa_order,
                                    enrichment_map) {
  
  res_list <- list()
  
  # 和原来一样，直接遍历向量
  for (target in intersect_merge) {
    
    virus_vec <- virus_exp_table[rownames(virus_exp_table) == target, ]
    bact_vec <- bact_exp_table[rownames(bact_exp_table) == target, ]
    
    if (length(virus_vec) == 0 | length(bact_vec) == 0) next
    
    df <- data.frame(
      Sample_ID = colnames(virus_exp_table),
      virus = as.numeric(virus_vec),
      bact = as.numeric(bact_vec)
    ) %>%
      left_join(sample_info %>% select(Sample_ID, Group_Area), by = "Sample_ID")
    
    for (i in seq_along(group_levels)) {
      
      stage_data <- df %>% filter(Group_Area == group_levels[i])
      
      if (nrow(stage_data) < 6) next
      
      ct <- cor.test(stage_data$virus, stage_data$bact, method = "spearman")
      
      res_list[[length(res_list) + 1]] <- data.frame(
        Stage = stage_labels[i],
        r = ct$estimate,
        p_value = ct$p.value,
        taxa = target
      )
    }
  }
  
  correlation_result <- bind_rows(res_list) %>%
    left_join(enrichment_map, by = "taxa") %>%   # 唯一新增：join 组别信息
    mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
    filter(p_value < 0.05) %>%
    mutate(
      neg_log_p = -log10(p_value),
      Stage = factor(Stage, levels = stage_labels),
      taxa = factor(taxa, levels = taxa_order),
      Enriched_Group = factor(Enriched_Group, levels = names(intersect_merge_list))
    )
  
  ggplot(correlation_result, aes(x = Stage, y = taxa, fill = r, size = neg_log_p)) +
    geom_point(shape = 21, alpha = 0.85, color = "black", stroke = 0.3) +
    scale_fill_gradient2(
      low = "#6e48fb", mid = "white", high = "#ffa500",
      midpoint = 0, limits = c(-1, 1), name = "Spearman r"
    ) +
    scale_size_continuous(range = c(1, 6), name = expression(-log[10](p))) +
    facet_grid(Enriched_Group ~ ., scales = "free_y", space = "free_y",labeller = labeller(
      Enriched_Group = function(x) paste0("Enriched in\n", x) ))  +
    labs(x = "Developmental Stage", y = "Taxa") +
    theme_classic(base_size = 8) +
    theme(
      strip.text.y = element_text(angle = 270, hjust=0.5, vjust = 0,face = "bold"),
      strip.background = element_rect(fill = "grey90", color = NA),
      panel.spacing.y = unit(0.5, "lines")
    )
}
# ==================== 5. 普氏分析函数 ====================
procrustes_analysis <- function(virus_filter, bact_filter, plot_title) {
  virus_dist <- vegdist(virus_filter, method = "bray")
  bact_dist <- vegdist(bact_filter, method = "bray")
  
  set.seed(123)
  mds_virus <- monoMDS(virus_dist)
  mds_bact <- monoMDS(bact_dist)
  
  pro_result <- procrustes(mds_virus, mds_bact, symmetric = TRUE)
  pro_test <- protest(mds_virus, mds_bact, permutations = 999)
  
  pro_points <- data.frame(
    bact_x = pro_result$Yrot[, 1],
    bact_y = pro_result$Yrot[, 2],
    virus_x = pro_result$X[, 1],
    virus_y = pro_result$X[, 2],
    sample = rownames(pro_result$Yrot)
  )
  
  pro_points_long <- rbind(
    pro_points %>% select(x = virus_x, y = virus_y, sample) %>% mutate(Type = "virus"),
    pro_points %>% select(x = bact_x, y = bact_y, sample) %>% mutate(Type = "bacteria")
  )
  
  p <- ggplot() +
    geom_segment(data = pro_points,
                 aes(x = virus_x, y = virus_y, xend = bact_x, yend = bact_y),
                 arrow = arrow(length = unit(0.2, "cm")),
                 color = "gray70", alpha = 0.8) +
    geom_point(data = pro_points_long,
               aes(x = x, y = y, color = Type, shape = Type),
               alpha = 0.7, size = 3) +
    scale_color_manual(values = c("virus" = "#E41A1C", "bacteria" = "#377EB8")) +
    scale_shape_manual(values = c("virus" = 17, "bacteria" = 19)) +
    labs(title = plot_title, x = "Dimension 1", y = "Dimension 2") +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(size = 7),
      axis.text = element_text(size = 6),
      plot.title = element_text(size = 8),
      legend.position = c(0.85, 0.15)
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
    annotate("text", x = -Inf, y = Inf,
             label = paste("M² =", round(pro_test$ss, 3), "p =", pro_test$signif),
             hjust = -0.05, vjust = 1.5, size = 3)
  
  return(p)
}

# ==================== 6. 多样性相关性函数（更新为使用 Group_Area）====================
diversity_correlation <- function(shannon_data, group_area, color) {
  # 从 Group_Area 中提取 Area 和 Group
  parts <- strsplit(group_area, "_")[[1]]
  area <- parts[1]
  group <- parts[2]
  
  bac <- shannon_data[shannon_data$Area == area & shannon_data$Group == group & shannon_data$type == "Bacteria", ]
  bac <- bac[, c("Sample_ID", "diversity")]
  names(bac)[2] <- "Bacteria"
  
  vir <- shannon_data[shannon_data$Area == area & shannon_data$Group == group & shannon_data$type == "Viruse", ]
  vir <- vir[, c("Sample_ID", "diversity")]
  names(vir)[2] <- "Viruse"
  
  merge_data <- merge(bac, vir, by = "Sample_ID")
  merge_data$Sample_ID <- NULL
  
  cor.result <- corr.test(merge_data, method = "pearson", adjust = "none")
  r <- cor.result$r[1, 2]
  p <- cor.result$p[1, 2]
  
  stat_text <- paste0("r = ", round(r, 4), ", p = ", round(p, 4))
  title_text <- paste(group_area, "Alpha Diversity Correlation")
  
  p_plot <- ggplot(merge_data, aes(x = Bacteria, y = Viruse)) +
    geom_point(shape = 16, size = 1, color = color) +
    stat_smooth(method = lm, se = FALSE, color = "black", size = 0.5) +
    theme_bw() +
    labs(title = title_text, x = "Bacteria", y = "Viruse") +
    theme(
      axis.text = element_text(size = 6),
      axis.title = element_text(size = 7),
      axis.line = element_line(color = "black"),
      panel.grid = element_blank(),
      plot.title = element_text(size = 8)
    ) +
    annotate("text", x = min(merge_data$Bacteria, na.rm = TRUE),
             y = max(merge_data$Viruse, na.rm = TRUE),
             label = stat_text, hjust = 0, vjust = 1, size = 2, color = "black")
  
  return(p_plot)
}

message("✓ Public utility functions loaded")