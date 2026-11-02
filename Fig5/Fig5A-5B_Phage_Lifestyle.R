# ============================================================
# Figure 5A-5B: 噬菌体生活方式与前噬菌体维持策略
# ============================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggsignif)
  library(dplyr)
  library(scales)
})

source("../R/00_config.R")

dir.create("../output/Fig5", recursive = TRUE, showWarnings = FALSE)

# ==================== 1. Figure 5A: 噬菌体生活方式比例 ====================
phastyle <- fread("../data/blast_viruse.tsv")
# phastyle$
setDT(phastyle )
temperate_prop <- phastyle[, .(temperate_prop = sum(predicted_label == "Temperate") / .N),
                           by = Sample]
temperate_prop <- temperate_prop %>%
  filter(temperate_prop > 0)
# temperate_prop[, temperate_log2 := log2(temperate_prop + 1)]

meta <- fread("../data/meta.tsv")
meta[, Group_Area := paste(Area, Group, sep = "_")]

df <- merge(temperate_prop, meta, by.x = "Sample", by.y = "Sample_ID")
df <- df[Group_Area %in% TARGET_GROUPS]

p_5a <- ggplot(df, aes(x = Group, y = temperate_prop, fill = Group)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.85, width = 0.7) +
  facet_grid(~ Area, scales = "free", space = "free") +
  scale_fill_manual(values = group_color) +
  labs(
    x = NULL,
    y = expression("Temperate Phage Proportion")
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "#F5F5F5", color = NA),
    plot.margin = margin(15, 10, 10, 10)
  ) +
  # 直接在这里添加显著性标记
  lapply(names(all_comparisons), function(area) {
    geom_signif(
      data = subset(df, Area == area),
      comparisons = all_comparisons[[area]],
      test = "wilcox.test",
      # 关键：只有 p < 0.05 才显示星号，否则显示空字符串（不标注）
      map_signif_level = function(p) {
        ifelse(p < 0.001, "***",
               ifelse(p < 0.01, "**",
                      ifelse(p < 0.05, "*", "")))
      },
      textsize = 3,
      step_increase = 0.12,
      tip_length = 0.01
    )
  })

ggsave("../output/Fig5/Figure5A_Lifestyle.pdf", p_5a, width = 4, height = 4)
message("✓ Figure 5A saved")

# ==================== 2. Figure 5B: 温和噬菌体比例对比 ====================
# 加载Fig4的中间结果
dt <- readRDS("../output/intermediate/shared_events_dt.rds")
final_table <- readRDS("../output/intermediate/final_table.rds")

anno <- fread("../data/blast_viruse4.tsv")
anno[, Type := ifelse(tolower(predicted_label) == "temperate", "Temperate", "Lytic")]
anno_phage <- anno[, .(Contig_name, Type, Genus_anno = host_annotation2)]

shared_base <- dt[, .(Contigs, Groups, Paired_Type, Genus = Annotation_cluster)]
shared_base <- shared_base[!is.na(Genus) & Genus != ""]
shared <- merge(shared_base, anno_phage, by.x = "Contigs", by.y = "Contig_name", all.x = TRUE)
shared <- shared[!is.na(Type)]

# 背景比例
bg_raw <- anno[, .(Genus = host_annotation2, Type)][!is.na(Genus) & Genus != "",
                                                    .N, by = .(Genus, Type)]
bg_wide <- dcast(bg_raw, Genus ~ Type, value.var = "N", fill = 0)
setDT(bg_wide)
bg_wide[, total := Lytic + Temperate]
bg_wide[, prop_bg := round(Temperate / total * 100, 1)]
bg <- bg_wide[, .(Genus, prop_bg, bg_total = total)]

# 共享比例
shared_raw <- shared[, .N, by = .(Genus, Groups, Type, Paired_Type)]
shared_wide <- dcast(shared_raw, Genus + Groups + Paired_Type ~ Type, value.var = "N", fill = 0)
setDT(shared_wide)
shared_wide[, total := Lytic + Temperate]
shared_wide[, prop_shared := round(Temperate / total * 100, 1)]
shared_prop <- shared_wide[, .(Genus, Groups, Paired_Type, prop_shared, n_shared = total)]
shared_prop <- shared_prop[n_shared >= 5]

# 合并并计算差异
comp <- merge(shared_prop, bg, by = "Genus", all.x = TRUE)
comp[, diff := prop_shared - fifelse(is.na(prop_bg), 50, prop_bg)]

# 使用核心属
key_genus <- unique(final_table$Genera[final_table$Venn_Label == "G"])
plot_df <- comp[Genus %in% key_genus][, prop_bg := fifelse(is.na(prop_bg), 50, prop_bg)]
plot_df[, genus_idx := .GRP, by = Genus]
set.seed(42)
plot_df[, jitter_x := genus_idx + (runif(.N, -0.25, 0.25))]

plot_df[, Transmission := fifelse(grepl("M0d|M21d", Groups),
                                  "Mother to offspring",
                                  "Oral to meconium")]
plot_df[, Transmission := factor(Transmission, levels = c("Mother to offspring", "Oral to meconium"))]

col_paired   <- "#C0392B"
col_unpaired <- "#2E6FA3"
col_bg       <- "#8C8C8C"

bg_ref <- unique(plot_df[, .(genus_idx, Genus, prop_bg)])
notable <- plot_df[abs(diff) >= 12]
notable[, vj := fifelse(diff > 0, -0.9, 1.7)]

p_5b <- ggplot(plot_df) +
  geom_segment(data = bg_ref,
               aes(x = genus_idx - 0.3, xend = genus_idx + 0.3, y = prop_bg, yend = prop_bg),
               color = col_bg, linewidth = 1.5, alpha = 0.45) +
  geom_point(data = bg_ref, aes(x = genus_idx, y = prop_bg),
             color = col_bg, size = 2.6, shape = 18, alpha = 0.9) +
  geom_segment(aes(x = jitter_x, xend = genus_idx, y = prop_shared, yend = prop_bg),
               color = "#C9C9C9", linewidth = 0.25, alpha = 0.5) +
  geom_point(aes(x = jitter_x, y = prop_shared, fill = Paired_Type,
                 shape = Transmission, size = n_shared),
             colour = "#333333", stroke = 0.5, alpha = 0.92) +
  geom_text(data = notable,
            aes(x = jitter_x, y = prop_shared,
                label = paste0(round(prop_shared, 0), "%"), vjust = vj),
            size = 2.8, fontface = "bold", color = "#2D2D2D", check_overlap = TRUE) +
  scale_fill_manual(values = c("Paired" = col_paired, "Unpaired" = col_unpaired),
                    name = "Pairing", drop = FALSE) +
  scale_shape_manual(values = c("Mother to offspring" = 21, "Oral to meconium" = 24),
                     name = "Transmission route", drop = FALSE) +
  scale_size_continuous(range = c(2, 5.5), name = "Contig count") +
  scale_x_continuous(breaks = unique(plot_df$genus_idx), labels = unique(plot_df$Genus),
                     expand = expansion(mult = 0.08)) +
  scale_y_continuous(labels = label_percent(scale = 1), limits = c(0, 100),
                     breaks = seq(0, 100, 25), expand = expansion(mult = c(0.02, 0.06))) +
  guides(
    shape = guide_legend(order = 1, override.aes = list(fill = "#666666", colour = "#333333", size = 3)),
    fill  = guide_legend(order = 2, override.aes = list(shape = 21, colour = "#333333")),
    size  = guide_legend(order = 3)
  ) +
  labs(x = NULL, y = "Temperate phage proportion (%)") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "#F0F0F0", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "#D9D9D9", fill = NA, linewidth = 0.5),
    axis.text.x = element_text(angle = 35, hjust = 1, size = 11, face = "italic"),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 11, face = "bold", margin = margin(r = 8)),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.box.just = "center"
  )

ggsave("../output/Fig5/Figure5B_Temperate_Proportion.pdf", p_5b, width = 5, height = 5)
message("✓ Figure 5B saved")
message("✓ Phage lifestyle analysis completed")





##############5c############
Integrase1 <- read.csv("D:/组会/进展汇报/FIG3_2/FIG3_2/for_Integrase_abundance_matrix.tsv", sep="\t", check.names = F)
Integrase2 <- read.csv("D:/组会/进展汇报/FIG3_2/FIG3_2/for_Integrase_abundance_matrix2.tsv", sep="\t", check.names = F)
Integrase <- merge(Integrase1, Integrase2, by="contig", all = T)

Integrase3 <- read.csv("D:/组会/进展汇报/FIG3_2/FIG3_2/for_recombinase_abundance_matrix.tsv", sep="\t", check.names = F)
Integrase4 <- read.csv("D:/组会/进展汇报/FIG3_2/FIG3_2/for_recombinase_abundance_matrix2.tsv", sep="\t", check.names = F)
Recombinase <- merge(Integrase3, Integrase4, by="contig", all = T)

all_data <- merge(Integrase, Recombinase, by="contig", all = T)
# all_data <-Recombinase
rownames(all_data) <- all_data$contig
all_data <- all_data[, colnames(all_data) %in% c("contig") == F]

simple_combine <- function(df1, df2) {
  rownames(df1) <- df1$contig
  rownames(df2) <- df2$contig
  df1 <- df1[, colnames(df1) != "contig", drop = FALSE]
  df2 <- df2[, colnames(df2) != "contig", drop = FALSE]
  
  duplicate_cols <- intersect(colnames(df1), colnames(df2))
  all_rows <- unique(c(rownames(df1), rownames(df2)))
  all_cols <- unique(c(colnames(df1), colnames(df2)))
  result <- matrix(0, nrow = length(all_rows), ncol = length(all_cols))
  rownames(result) <- all_rows
  colnames(result) <- all_cols
  result[rownames(df1), colnames(df1)] <- as.matrix(df1)
  result[rownames(df2), colnames(df2)] <- result[rownames(df2), colnames(df2)] + as.matrix(df2)
  
  return(result)
}
Integrase_combined <- simple_combine(Integrase1, Integrase2)
Recombinase_combined <- simple_combine(Integrase3, Integrase4)
all_data <- simple_combine(
  data.frame(contig = rownames(Integrase_combined), Integrase_combined, check.names = FALSE),
  data.frame(contig = rownames(Recombinase_combined), Recombinase_combined, check.names = FALSE)
)
abundance <- function(mat, min_sum = 20) {
  mat <- mat[, order(colnames(mat))]
  mat[is.na(mat)] <- 0
  keep_rows <- rowSums(mat) >= min_sum
  mat <- mat[keep_rows, , drop = FALSE]
  apply(mat, 2, function(x) x / sum(x))
}
all_data[is.na(all_data)] <- 0
all_data <- abundance(all_data)

library(dplyr)
library(tidyr)
library(tidyverse)
blast_viruse <- read.csv("../data/blast_viruse4.tsv", header=TRUE,sep="\t")
result <- data.frame(all_data, check.names = F) %>%
  rownames_to_column("Integrase_id") %>%
  pivot_longer(-Integrase_id, names_to = "sample", values_to = "abundance") %>%
  filter(abundance > 0) %>%
  mutate(contig_id = str_remove(Integrase_id, "_[^_]*$")) %>%
  left_join(blast_viruse, by = c("contig_id" = "Contig_name")) %>%
  mutate(genus = host_annotation2) %>%
  group_by(genus, sample) %>%
  summarise(total_abundance = sum(abundance, na.rm = TRUE), .groups = "drop")
result2 <- data.frame(all_data, check.names = F) %>%
  rownames_to_column("Integrase_id") %>%
  pivot_longer(-Integrase_id, names_to = "sample", values_to = "abundance") %>%
  filter(abundance > 0) %>%
  mutate(contig_id = str_remove(Integrase_id, "_[^_]*$")) %>%
  left_join(blast_viruse, by = c("contig_id" = "Contig_name")) %>%
  mutate(genus = host_annotation2) %>%
  group_by(genus, sample)
result2 <- merge(result2[c("abundance","sample","genus")], meta, by.x = "sample", by.y = "Sample_ID")
# result2 <- result2[result2$Group=="CME",]
result <- merge(result, meta, by.x = "sample", by.y = "Sample_ID")

target_genus <-c(  "Corynebacterium", "Glutamicibacter", "Microbacterium","Psychrobacter",  "Aerococcus")

genus <-target_genus
# genus <-unique(final_table$Genera[final_table$Venn_Label=="G"])
result <- result[result$genus %in% genus, ]

library(ggplot2)
library(ggsignif)
library(dplyr)
result <- result[result$Group_Area%in%GROUP_LEVELS,]
p <- ggplot(result, aes(x = Group_Area, y = total_abundance, fill = Group_Area)) +
  geom_boxplot() +
  facet_wrap(~ genus, scales = "free",nrow=1) +
  theme_bw() +
  scale_fill_manual(values = group_color) +
  lapply(unique(result$genus), function(g) {
    current_data <- result %>% filter(genus == g)
    present_groups <- unique(current_data$Group_Area)
    valid_comps <- list()
    for (pair in comparisons) {
      if (all(pair %in% present_groups)) {
        valid_comps[[length(valid_comps) + 1]] <- pair
      }
    }
    if (length(valid_comps) > 0) {
      y_max <- max(current_data$total_abundance, na.rm = TRUE)
      y_min <- min(current_data$total_abundance, na.rm = TRUE)
      lapply(seq_along(valid_comps), function(i) {
        geom_signif(
          data = current_data,
          comparisons = list(valid_comps[[i]]),
          test = "wilcox.test",
          map_signif_level = function(p) {
            ifelse(p < 0.001, "***",
                   ifelse(p < 0.01, "**",
                          ifelse(p < 0.05, "*", NA)))
          },
          tip_length = 0.01,
          textsize = 4,
          vjust = 0.5,
          y_position = y_max + (y_max - y_min) * 0.1 * i
        )
      })
    }
  })+
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "#F5F5F5", color = NA),
    plot.margin = margin(15, 10, 10, 10)
  )
ggsave("../output/Fig5/Figure5C_integase.pdf", p, width =10, height = 5)
print(p)
library(tidyverse)
library(ggpubr)
desc_table <- result %>%
  group_by(genus, Group_Area) %>%
  summarise(
    n = n(),
    mean = mean(total_abundance, na.rm = TRUE),
    sd = sd(total_abundance, na.rm = TRUE),
    median = median(total_abundance, na.rm = TRUE),
    .groups = "drop"
  )
sig_table <- result %>%
  group_by(genus) %>%
  rstatix::wilcox_test(total_abundance ~ Group_Area) %>%
  rstatix::adjust_pvalue(method = "BH") %>%
  rstatix::add_significance() %>%
  select(genus, group1, group2, p.adj, p.adj.signif)