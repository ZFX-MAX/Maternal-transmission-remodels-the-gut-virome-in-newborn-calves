# ============================================================
# Figure 5C-5D: 前噬菌体基因图谱与测序深度分析
# ============================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(gggenes)
  library(patchwork)
  library(zoo)
  library(fs)
})

source("../R/00_config.R")

dir.create("../output/Fig5", recursive = TRUE, showWarnings = FALSE)

BASE_DIR <- "../data/prophage_detection"
DEPTH_FILE <- file.path(BASE_DIR, "ALL_DNA_depth.tsv")

# ==================== 1. 模块注释规则 ====================
MODULE_RULES <- list(
  "Integration"    = c("integrase", "recombinase", "site-specific"),
  "Packaging"      = c("terminase", "packaging", "portal"),
  "Head"           = c("head", "capsid", "scaffold", "core", "coat"),
  "Tail"           = c("tail", "fiber", "baseplate", "spike", "tube"),
  "Lysis"          = c("lysin", "holin", "endolysin", "muramidase", "peptidoglycan"),
  "DNA Metabolism" = c("polymerase", "exonuclease", "helicase", "primase", "ligase"),
  "Regulation"     = c("repressor", "regulator", "Cro", "CI", "antiterminator")
)

assign_module <- function(desc, rules_list) {
  desc_lower <- tolower(desc)
  for(module_name in names(rules_list)) {
    if(any(str_detect(desc_lower, regex(paste(rules_list[[module_name]], collapse = "|"), ignore_case = TRUE)))) {
      return(module_name)
    }
  }
  return("Hypothetical")
}

# ==================== 2. GFF解析与注释加载 ====================
load_genus_annotation <- function(genus, base_dir) {
  diamond_file <- file.path(base_dir, genus, sprintf("%s_bact_prophage_nr.tsv", genus))
  nr_file <- file.path(base_dir, genus, "nr_desc.tsv")
  
  if(!file.exists(diamond_file)) return(tibble(gene_id = character(), desc = character()))
  
  diamond <- read_tsv(diamond_file,
                      col_names = c("gene_id", "sseqid", "evalue", "bitscore"),
                      show_col_types = FALSE)
  
  nr_desc <- tibble(sseqid = character(), description = character())
  if(file.exists(nr_file) && file.info(nr_file)$size > 0) {
    lines <- read_lines(nr_file)
    lines <- lines[!grepl("^#", lines)]
    lines <- lines[trimws(lines) != ""]
    
    if(length(lines) > 0) {
      nr_desc <- tibble(line = lines) %>%
        mutate(
          sseqid = str_extract(line, "^[^\\s]+"),
          description = str_replace(line, "^[^\\s]+\\s+", "") %>% str_trim()
        ) %>%
        select(sseqid, description) %>%
        distinct(sseqid, .keep_all = TRUE)
    }
  }
  
  result <- diamond %>%
    left_join(nr_desc, by = "sseqid") %>%
    mutate(
      is_hypothetical = case_when(
        is.na(description) ~ TRUE,
        str_detect(tolower(description), "hypothetical|unknown|uncharacterized|putative") ~ TRUE,
        TRUE ~ FALSE
      ),
      desc = if_else(!is.na(description) & !is_hypothetical,
                     str_replace(description, "\\s*\\[.*\\]\\s*$", "") %>% str_trim(),
                     "hypothetical protein")
    ) %>%
    arrange(gene_id, is_hypothetical, desc(bitscore), evalue) %>%
    group_by(gene_id) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    select(gene_id, desc)
  
  return(result)
}

parse_gff_genes <- function(gff_path, genus_name, base_dir) {
  if(!file.exists(gff_path)) return(tibble())
  
  genes <- read_tsv(gff_path, comment = "#", col_names = FALSE, show_col_types = FALSE) %>%
    filter(X3 == "CDS") %>%
    mutate(
      Genus   = genus_name,
      contig  = X1,
      start   = as.integer(X4),
      end     = as.integer(X5),
      strand  = ifelse(X7 == "+", 1L, -1L),
      num     = str_extract(X9, "ID=(?:.*_)?(\\d+)", group = 1),
      gene_id = paste0(X1, "_", num),
      product = str_extract(X9, "product=([^;]+)", group = 1)
    ) %>%
    select(Genus, contig, start, end, strand, gene_id, product)
  
  annot <- load_genus_annotation(genus_name, base_dir)
  
  if(nrow(genes) > 0 && nrow(annot) > 0) {
    genes %>%
      left_join(annot, by = "gene_id") %>%
      mutate(
        desc = coalesce(desc, product, gene_id, "hypothetical protein"),
        desc = str_trim(desc),
        module = map_chr(desc, ~assign_module(.x, MODULE_RULES))
      )
  } else if(nrow(genes) > 0) {
    genes %>%
      mutate(
        desc = coalesce(product, gene_id, "hypothetical protein"),
        module = "Hypothetical"
      )
  } else {
    genes
  }
}

# 加载所有基因
gff_files <- dir_ls(BASE_DIR, recurse = TRUE, glob = "*_bact_prophage.gff")
genus_names <- path_file(path_dir(gff_files))
GENUS_CONFIG <- setNames(as.list(path_rel(gff_files, start = ".")), genus_names)

message("Parsing GFF files...")
ALL_GENES <- bind_rows(
  lapply(names(GENUS_CONFIG), function(genus) {
    parse_gff_genes(GENUS_CONFIG[[genus]], genus, base_dir = BASE_DIR)
  })
) %>%
  mutate(
    gene_len = end - start + 1,
    module = str_trim(if_else(is.na(module) | module == "", "Hypothetical", module))
  ) %>%
  arrange(Genus, contig, start)

# ==================== 3. 前噬菌体图谱绘制函数 ====================
plot_prophage_combined <- function(contig_name, prophage_start, prophage_end,
                                   attL, attR, genes_df, depth_file, meta_df,
                                   save_path = NULL, width = 10, height = 6) {
  genes_target <- genes_df %>% filter(contig == contig_name) %>% arrange(start)
  seq_len <- if(nrow(genes_target) > 0) max(genes_target$end) + 2000 else prophage_end + 5000
  x_limits <- c(0, seq_len)
  
  mod_colors <- c(
    "Integration" = "#8e44ad", "Packaging" = "#c0392b",
    "Head" = "#2980b9", "Tail" = "#27ae60",
    "Lysis" = "#d35400", "DNA Metabolism" = "#16a085",
    "Regulation" = "#2c3e50", "Hypothetical" = "#bdc3c7"
  )
  
  # 基因图谱
  p_genes <- ggplot(genes_target) +
    annotate("rect", xmin = prophage_start, xmax = prophage_end,
             ymin = 0.5, ymax = 1.5, fill = "gray85", alpha = 0.4) +
    geom_gene_arrow(
      aes(xmin = start, xmax = end, y = 1, fill = module, forward = strand == 1),
      arrowhead_height = unit(4, "mm"),
      arrow_body_height = unit(3, "mm"),
      size = 0.3, color = "black", alpha = 0.95
    ) +
    annotate("segment", x = attL, xend = attL, y = 0.5, yend = 1.5,
             color = "#e74c3c", linetype = "dashed", size = 0.6) +
    annotate("segment", x = attR, xend = attR, y = 0.5, yend = 1.5,
             color = "#e74c3c", linetype = "dashed", size = 0.6) +
    scale_fill_manual(values = mod_colors, na.value = "#bdc3c7",
                      name = "Module", drop = FALSE) +
    scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.5, 1.5), expand = c(0, 0)) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      axis.line.x.bottom = element_line(color = "gray", linewidth = 0.8, linetype = "dotted"),
      legend.position = "right",
      plot.margin = margin(10, 45, 5, 10, "mm")
    )
  
  # 测序深度
  depth <- read_tsv(depth_file, show_col_types = FALSE)
  colnames(depth)[colnames(depth) == "Pos"] <- "Position"
  d <- depth %>% filter(Contig == contig_name) %>% select(Sample, Position, Depth)
  
  if(nrow(d) > 0) {
    d_grouped <- d %>%
      mutate(Window_Mid = ((Position - 1) %/% 500) * 500 + 250) %>%
      left_join(meta_df, by = c("Sample" = "Sample_ID")) %>%
      filter(Group_Area %in% TARGET_GROUPS) %>%
      group_by(Area, Group_Area, Window_Mid) %>%
      summarise(Mean_Depth = mean(Depth, na.rm = TRUE),
                SD_Depth = sd(Depth, na.rm = TRUE),
                N = n(), .groups = "drop") %>%
      filter(N >= 3) %>%
      group_by(Area, Group_Area) %>%
      mutate(Mean_Depth = rollapply(Mean_Depth, width = 5, FUN = mean,
                                    fill = NA, align = "center", na.rm = TRUE)) %>%
      ungroup() %>%
      filter(!is.na(Mean_Depth), Mean_Depth > 0)
    
    p_depth <- ggplot(d_grouped, aes(x = Window_Mid, y = Mean_Depth,
                                     color = Group_Area, fill = Group_Area)) +
      geom_line(linewidth = 1.2) +
      geom_ribbon(aes(ymin = pmax(0.1, Mean_Depth - SD_Depth),
                      ymax = Mean_Depth + SD_Depth), alpha = 0.12, color = NA) +
      scale_color_manual(values = group_color, guide = "none") +
      scale_fill_manual(values = group_color, guide = "none") +
      scale_y_continuous(trans = "log10", name = "Depth (log10)",
                         breaks = c(1, 10, 50, 100, 500),
                         labels = c("1", "10", "50", "100", "500")) +
      coord_cartesian(ylim = c(0.5, NA)) +
      scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
      facet_wrap(~Area, strip.position = "right") +
      theme_minimal(base_size = 10) +
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background = element_rect(fill = "gray95", color = NA),
        strip.text = element_text(face = "bold", size = 11, angle = 0),
        strip.placement = "outside",
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = margin(0, 45, 5, 10, "mm")
      )
  } else {
    p_depth <- ggplot() + theme_void()
  }
  
  combined <- p_genes / p_depth + plot_layout(heights = c(1, 2.5))
  
  if(!is.null(save_path)) {
    ggsave(save_path, combined, width = width, height = height, dpi = 300, bg = "white")
  }
  
  return(combined)
}

# ==================== 4. 绘制代表性前噬菌体 ====================
meta <- fread("../data/meta.tsv")
meta[, Group_Area := paste(Area, Group, sep = "_")]

# Fig5C: Aerococcus
p_5c <- plot_prophage_combined(
  contig_name = "0dKG00334.NODE_1132_length_21712_cov_80.634114",
  prophage_start = 8118, prophage_end = 18862,
  attL = 8118, attR = 18862,
  genes_df = ALL_GENES, depth_file = DEPTH_FILE, meta_df = meta,
  save_path = "../output/Fig5/Figure5C_Aerococcus_Prophage.pdf"
)
message("✓ Figure 5C saved")

# Fig5D: Psychrobacter
p_5d <- plot_prophage_combined(
  contig_name = "0dKG00188.NODE_150_length_29412_cov_48.624451",
  prophage_start = 23000, prophage_end = 27000,
  attL = 23000, attR = 27000,
  genes_df = ALL_GENES, depth_file = DEPTH_FILE, meta_df = meta,
  save_path = "../output/Fig5/Figure5D_Psychrobacter_Prophage.pdf"
)
message("✓ Figure 5D saved")
message("✓ Prophage mapping completed")

