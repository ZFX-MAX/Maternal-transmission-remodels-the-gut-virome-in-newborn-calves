# ============================================================
# Figure 6C: AMG Pathway Analysis and Metabolic Schematic
# Description: Analyze phage-encoded AMGs, compare with bacterial
#              pathways, and generate metabolic schematic diagrams
# Dependencies: data.table, tidyverse, ggplot2, patchwork, stringr
# ============================================================
suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(stringr)
})

source("../R/00_config.R")
dir.create("../output/Fig6", recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. Load KEGG Mapping Files (格式已根据真实文件修正)
# ============================================================
message("Loading KEGG mapping files...")

# ko_name.txt: 第一列是 K00001 (无 ko: 前缀), 第二列是描述
ko_name <- fread("../data/kegg/ko_name.txt", header = FALSE,
                 col.names = c("KO", "KO_description"), sep = "\t")
# 不需要 gsub，因为本来就没有 ko: 前缀

# pathway_name.txt: 第一列是 map01100 (map 前缀), 第二列是描述
pathway_name <- fread("../data/kegg/pathway_name.txt", header = FALSE, sep = "\t",
                      col.names = c("Pathway_number", "Pathway_description"))
pathway_name[, Pathway_number := gsub("^map", "ko", Pathway_number)]  # map -> ko

# ko_pathway_name.txt: 两列都有前缀 (ko:K00001, path:map00010)
ko2pathway <- fread("../data/kegg/ko_pathway_name.txt", header = FALSE,
                    col.names = c("KO_raw", "Pathway_raw"))
ko2pathway[, KO := gsub("^ko:", "", KO_raw)]
ko2pathway[, Pathway := gsub("^path:map", "ko", Pathway_raw)]  # path:map -> ko
ko2pathway[, c("KO_raw", "Pathway_raw") := NULL]
ko2pathway <- unique(ko2pathway[!grepl("^ko01", Pathway)])  # 排除全局通路
ko2pathway <- merge(ko2pathway, pathway_name, by.x = "Pathway", by.y = "Pathway_number")
ko2pathway <- merge(ko2pathway, ko_name, by = "KO")

message("✓ KEGG mapping files loaded")

# ============================================================
# 2. Genus-Level Pathway Preference Analysis
# ============================================================
analyze_genus_preference <- function(df, core_genera, meta, target_groups, dataset_name) {
  desc_col <- intersect(c("Pathway_description", "KO_description"), colnames(df))
  desc_col <- if (length(desc_col) > 0) desc_col[1] else NULL
  
  id_cols <- c("Host_Genus", "Pathway_name")
  if (!is.null(desc_col)) id_cols <- c(id_cols, desc_col)
  
  all_sample_cols <- setdiff(colnames(df), id_cols)
  target_samples <- meta$Sample_ID[meta$Group %in% target_groups]
  sample_cols <- intersect(all_sample_cols, target_samples)
  
  message(sprintf("[%s] 筛选样本: %d 个 (目标组: %s)", 
                  dataset_name, length(sample_cols), paste(target_groups, collapse = "+")))
  
  dt <- melt.data.table(df, id.vars = id_cols, measure.vars = sample_cols,
                        variable.name = "Sample", value.name = "Abundance")
  
  dt[, Genus_Total := sum(Abundance, na.rm = TRUE), by = c("Sample", "Host_Genus")]
  dt[, Genus_Relative := Abundance / Genus_Total]
  dt[, Sample_Total := sum(Abundance, na.rm = TRUE), by = "Sample"]
  
  ko_by <- c("Sample", "Pathway_name")
  if (!is.null(desc_col)) ko_by <- c(ko_by, desc_col)
  dt[, KO_Total := sum(Abundance, na.rm = TRUE), by = ko_by]
  dt[, Global_Relative := KO_Total / Sample_Total]
  
  group_by <- c("Host_Genus", "Pathway_name")
  if (!is.null(desc_col)) group_by <- c(group_by, desc_col)
  
  results <- dt[Host_Genus %in% core_genera, .(
    mean_genus = mean(Genus_Relative, na.rm = TRUE),
    mean_global = mean(Global_Relative, na.rm = TRUE),
    p_value = tryCatch(
      wilcox.test(Genus_Relative, Global_Relative, 
                  alternative = "greater", exact = FALSE)$p.value,
      error = function(e) NA
    ),
    n_samples = .N
  ), by = group_by]
  
  results[, fold_change := mean_genus / (mean_global + 1e-10)]
  results[, p_adj := p.adjust(p_value, method = "fdr"), by = Host_Genus]
  results[, Dataset := dataset_name]
  
  return(results)
}

core_genera <- c("Psychrobacter", "Corynebacterium", "Glutamicibacter",
                 "Aerococcus", "Microbacterium", "Nocardioides",
                 "Arthrobacter", "Pseudomonas", "Moraxella", "Brachybacterium")
target_groups <- c("C0dOR", "C0dME")

meta <- fread("../data/meta.tsv")

# ============================================================
# 3. Run Pathway Preference Analysis
# ============================================================
message("Running genus-level pathway preference analysis...")

# 病毒 AMG 通路
df_virus_pw <- fread("../data/kegg/AMG_Pathway_abundance_by_Genus.tsv")
colnames(df_virus_pw)[1:2] <- c("Host_Genus", "Pathway_name")
df_virus_pw <- merge(df_virus_pw, unique(pathway_name), by.x = "Pathway_name", by.y = "Pathway_number")

results_virus <- analyze_genus_preference(df_virus_pw, core_genera, meta, target_groups, "Virus_AMG")
sig_virus <- results_virus[p_adj < 0.05 & fold_change > 2]


library(ggplot2)
library(dplyr)

genus_levels <- sort(unique(sig_virus$Host_Genus))
pathway_levels <- sort(unique(sig_virus$Pathway_description))

genus_df <- data.frame(
  Host_Genus = factor(genus_levels, levels = genus_levels),
  y_genus = seq_along(genus_levels)
)

pathway_df <- data.frame(
  Pathway_description = factor(pathway_levels, levels = pathway_levels),
  y_pathway = seq_along(pathway_levels)
)

plot_edges <- sig_virus %>%
  mutate(p_adj_safe = ifelse(p_adj == 0, min(p_adj[p_adj > 0], na.rm = TRUE) / 10, p_adj),
         Host_Genus = factor(Host_Genus, levels = genus_levels)) %>%
  left_join(genus_df, by = "Host_Genus") %>%
  left_join(pathway_df, by = "Pathway_description")

nodes_left <- data.frame(
  x = 0, y = genus_df$y_genus, 
  label = genus_df$Host_Genus, 
  Host_Genus = genus_df$Host_Genus
)

nodes_right <- data.frame(
  x = 1, y = pathway_df$y_pathway, 
  label = pathway_df$Pathway_description
)
fig<- ggplot() +
  geom_segment(
    data = plot_edges,
    aes(x = 0, xend = 1, y = y_genus, yend = y_pathway, 
        color = Host_Genus, size = -log10(p_adj_safe)),
    alpha = 0.7, lineend = "round"
  ) +
  geom_point(data = nodes_left, aes(x = x, y = y, color = Host_Genus), size = 5) +
  geom_text(data = nodes_left, aes(x = x - 0.03, y = y, label = label), 
            hjust = 1, vjust = 0.5, size = 4, fontface = "italic") +
  geom_point(data = nodes_right, aes(x = x, y = y), size = 4, color = "gray") +
  geom_text(data = nodes_right, aes(x = x + 0.03, y = y, label = label), 
            hjust = 0, vjust = 0.5, size = 3.5) +
  scale_color_brewer(palette = "Set3", name = "Host Genus", 
                     guide = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  scale_size_continuous(range = c(0.5, 2.5), name = "-log10(p_adj)") +
  
  # 调整 X 轴范围，给左右两侧的长文本留出空间
  scale_x_continuous(limits = c(-0.6, 2.5)) +
  theme_void() +
  theme(
    legend.position = "right", # 图例放在右侧，避免挤占底部空间
    legend.box = "vertical",
    # 增加右侧边距，防止右侧过长的代谢通路名称被画布裁剪
    # plot.margin = margin(20, 150, 20, 20), 
    # plot.title = element_text(hjust = 0.5, size = 16, face = "bold", margin = margin(b = 20)),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11)
  ) +
  coord_cartesian(clip = "off") 
ggsave("../output/Fig6/FigS4_link.pdf", fig, width = 12, height = 8)

message("✓ Genus-level pathway preference analysis completed")

# ============================================================
# 4. Build Supplementary Table
# ============================================================
message("Building supplementary table...")

AMG_all <- fread("../data/kegg/AMG_anno2.txt")
AMG_oral <- AMG_all[!is.na(host_annotation2) & Type == "Oral"]

panel_genus <- list(
  ko00250 = "Corynebacterium",
  ko00790 = c("Psychrobacter", "Microbacterium"),
  ko04122 = "Psychrobacter"
)

rows <- rbindlist(lapply(names(panel_genus), function(pw) {
  gens <- panel_genus[[pw]]
  kos <- ko2pathway[Pathway == pw, KO]
  
  v <- AMG_oral[Pathway_name == pw & host_annotation2 %in% gens & AMG.KO %in% kos,
                .(Pathway = pw, KO = AMG.KO, Source = "Virus", Host_Genus = host_annotation2,
                  Lifestyle = paste(unique(predicted_label), collapse = "/"),
                  Group = paste(unique(Group), collapse = "/"), Niche = "Oral")]
  v <- unique(v)
  
  b <- df_bact_ko[Pathway_name %in% kos & Host_Genus %in% gens,
                  .(Pathway = pw, KO = Pathway_name, Source = "Bacteria", Host_Genus = Host_Genus)]
  b[, `:=`(Lifestyle = NA_character_, Group = NA_character_, Niche = NA_character_)]
  
  rbind(v, b, fill = TRUE)
}))

bact_key <- unique(rows[Source == "Bacteria", .(k = paste(Pathway, Host_Genus, KO))])
rows[Source == "Virus", Shared_with_host := 
       ifelse(paste(Pathway, Host_Genus, KO) %in% bact_key$k, "Yes", "No")]
rows[Source == "Bacteria", Shared_with_host := "—"]

rows <- merge(rows, unique(ko_name[, .(KO, KO_description)]), by = "KO", all.x = TRUE)
rows <- merge(rows, unique(pathway_name[, .(Pathway_number, Pathway_description)]),
              by.x = "Pathway", by.y = "Pathway_number", all.x = TRUE)
rows[, Source := factor(Source, levels = c("Virus", "Bacteria"))]

supp <- rows[, .(Pathway, Pathway_description, KO, KO_description, Source,
                 Shared_with_host, Host_Genus, Group, Niche, Lifestyle)]
setorder(supp, Pathway, Source, Host_Genus, KO)

# ============================================================
# 5. Fetch KEGG Reaction Definitions (WITH CACHE)
# ============================================================
message("Processing reaction definitions...")

# ko_rm.txt: 两列都有前缀 (ko:K00001, rn:R00623)
ko2rn <- fread("../data/kegg/ko_rm.txt", header = FALSE,
               col.names = c("ko_raw", "rn_raw"), sep = "\t")
ko2rn[, ko := gsub("^ko:", "", ko_raw)]
ko2rn[, rn := gsub("^rn:", "", rn_raw)]
ko2rn[, c("ko_raw", "rn_raw") := NULL]

target_kos <- unique(supp$KO)
target_rns <- unique(ko2rn[ko %in% target_kos, rn])

# ★ 缓存文件路径
CACHE_FILE <- "../output/Fig6/rxn_defs_cache.tsv"

if (file.exists(CACHE_FILE)) {
  message(sprintf("✓ Found cached reaction definitions: %s", CACHE_FILE))
  rxn_defs <- fread(CACHE_FILE)
  missing_rns <- setdiff(target_rns, rxn_defs$rn)
  
  if (length(missing_rns) > 0) {
    message(sprintf("  Found %d missing reactions, downloading...", length(missing_rns)))
    get_definition <- function(rn) {
      url <- paste0("https://rest.kegg.jp/get/", rn)
      txt <- tryCatch(readLines(url, warn = FALSE, encoding = "UTF-8"), error = function(e) NULL)
      if (is.null(txt)) return(data.table(rn = rn, definition = NA_character_))
      idx <- grep("^DEFINITION", txt)
      if (length(idx) == 0) return(data.table(rn = rn, definition = NA_character_))
      def_line <- sub("^DEFINITION\\s*", "", txt[idx[1]])
      k <- idx[1] + 1
      while (k <= length(txt) && grepl("^\\s{5,}", txt[k])) {
        def_line <- paste(def_line, trimws(txt[k]))
        k <- k + 1
      }
      data.table(rn = rn, definition = def_line)
    }
    new_defs <- rbindlist(lapply(missing_rns, function(rn) {
      Sys.sleep(0.25)
      get_definition(rn)
    }))
    rxn_defs <- rbind(rxn_defs, new_defs)
    fwrite(rxn_defs, CACHE_FILE, sep = "\t")
  } else {
    message("✓ All required reactions found in cache")
  }
} else {
  message("No cache found. Downloading reaction definitions from KEGG...")
  get_definition <- function(rn) {
    url <- paste0("https://rest.kegg.jp/get/", rn)
    txt <- tryCatch(readLines(url, warn = FALSE, encoding = "UTF-8"), error = function(e) NULL)
    if (is.null(txt)) return(data.table(rn = rn, definition = NA_character_))
    idx <- grep("^DEFINITION", txt)
    if (length(idx) == 0) return(data.table(rn = rn, definition = NA_character_))
    def_line <- sub("^DEFINITION\\s*", "", txt[idx[1]])
    k <- idx[1] + 1
    while (k <= length(txt) && grepl("^\\s{5,}", txt[k])) {
      def_line <- paste(def_line, trimws(txt[k]))
      k <- k + 1
    }
    data.table(rn = rn, definition = def_line)
  }
  rxn_defs <- rbindlist(lapply(target_rns, function(rn) {
    Sys.sleep(0.25)
    get_definition(rn)
  }))
  fwrite(rxn_defs, CACHE_FILE, sep = "\t")
  message(sprintf("✓ Downloaded %d reactions and saved to cache", nrow(rxn_defs)))
}

# 解析反应定义
split_def <- function(def) {
  if (is.na(def)) return(list(substrates = NA_character_, products = NA_character_, direction = NA_character_))
  if (grepl(" <=> ", def, fixed = TRUE)) {
    parts <- strsplit(def, " <=> ", fixed = TRUE)[[1]]
    dir <- "Reversible (<=>)"
  } else if (grepl(" => ", def, fixed = TRUE)) {
    parts <- strsplit(def, " => ", fixed = TRUE)[[1]]
    dir <- "Irreversible (=>)"
  } else if (grepl(" = ", def, fixed = TRUE)) {
    parts <- strsplit(def, " = ", fixed = TRUE)[[1]]
    dir <- "Unknown (=)"
  } else {
    return(list(substrates = def, products = NA_character_, direction = "Unknown"))
  }
  list(substrates = trimws(parts[1]),
       products = ifelse(length(parts) > 1, trimws(parts[2]), NA_character_),
       direction = dir)
}

parsed <- rbindlist(lapply(rxn_defs$definition, split_def))
rxn_defs <- cbind(rxn_defs, parsed)

ko_rxn_info <- merge(ko2rn[ko %in% target_kos], rxn_defs, by = "rn", all.x = TRUE)
ko_summary <- ko_rxn_info[, .(
  Reactions = paste(unique(na.omit(rn)), collapse = "; "),
  Substrates = paste(unique(na.omit(substrates)), collapse = " | "),
  Products = paste(unique(na.omit(products)), collapse = " | "),
  Direction = paste(unique(na.omit(direction)), collapse = "; ")
), by = ko]
setnames(ko_summary, "ko", "KO")

supp_final <- merge(supp, ko_summary, by = "KO", all.x = TRUE)
setcolorder(supp_final, c("Pathway", "Pathway_description", "KO", "KO_description",
                          "Substrates", "Products", "Direction", "Reactions",
                          "Source", "Shared_with_host", "Host_Genus", "Group", "Niche", "Lifestyle"))

fwrite(supp_final, "../output/Fig6/Supplementary_Table_AMG_Final.csv")
message("✓ Supplementary table saved")

# ============================================================
# 6. Generate Metabolic Schematic
# ============================================================
message("Generating metabolic schematic diagram...")

target_pw <- c("ko00250", "ko00790", "ko04122")
plot_dt <- supp_final[Pathway %in% target_pw]

plot_dt[, Gene := str_extract(KO_description, "^[^,;]+")]
plot_dt[, Gene := str_replace(Gene, ".,\\s", "")]
plot_dt[, Gene := str_extract(Gene, "^[^;]+")]

blacklist <- c("ATP", "ADP", "GTP", "GDP", "NAD\\+", "NADH", "NADP\\+", "NADPH",
               "H2O", "Oxygen", "Orthophosphate", "Diphosphate", "Triphosphate",
               "Formate", "Ammonia", "CO2", "H\\+", "Reduced acceptor", "Acceptor",
               "S-Adenosyl-L-methionine", "5'-Deoxyadenosine", "L-Methionine", "Glycolaldehyde")
blacklist_pattern <- paste0("\\b(", paste(blacklist, collapse = "|"), ")\\b")

clean_metabolites <- function(txt) {
  if (is.na(txt) || txt == "") return(NA_character_)
  first_rxn <- str_split(txt, "\\| ")[[1]][1]
  cleaned <- str_remove_all(first_rxn, blacklist_pattern)
  cleaned <- str_replace_all(cleaned, "\\s*\\+\\s*", " + ")
  cleaned <- str_trim(cleaned)
  cleaned <- str_replace(cleaned, "^\\+\\s*|\\s*\\+$", "")
  return(cleaned)
}

plot_dt[, Core_Substrates := sapply(Substrates, clean_metabolites)]
plot_dt[, Core_Products := sapply(Products, clean_metabolites)]
plot_dt[, is_viral_AMG := (Source == "Virus")]

plot_edges <- plot_dt[, .(
  Pathway, Host_Genus, KO, Gene,
  From = Core_Substrates, To = Core_Products,
  Source, Shared_with_host, is_viral_AMG
)]
plot_edges <- plot_edges[!is.na(From) & From != "" & !is.na(To) & To != ""]

fwrite(plot_edges, "../output/Fig6/Plot_Edges_Clean.csv")

# ============================================================
# 7. Draw Schematic Diagram
# ============================================================
N <- rbind(
  data.table(panel = 1, name = c("Oxaloacetate", "Aspartate", "Adenylosuccinate", "AMP", "DNA/RNA",
                                 "Asparagine", "Pyruvate", "Alanine"),
             x = c(0, 0, 0, 0, 0, 1.6, -1.6, -1.6), 
             y = c(4, 3, 2, 1, 0, 3, 4, 3),
             branch = c(F, TRUE, F, F, F, F, F, F)),
  data.table(panel = 2, name = c("GTP", "DHN-P3", "HMDP", "HMDP-PP", "Dihydropteroate", "DHF", "THF", "dTMP/purines"),
             x = 0, y = 7:0, 
             branch = c(F, F, F, F, F, TRUE, F, F)),
  data.table(panel = 3, name = c("Cysteine", "S-persulfide", "CysO-S-cysteine", "L-Cysteine",
                                 "Fe-S/other S sinks", "GTP/cPMo", "Precursor Z", "Molybdopterin", "Moco", "NarG/FDH"),
             x = c(0, 0, 0, 0, -1.8, 2, 2, 2, 2, 2), 
             y = c(4, 3, 2, 1, 3, 4, 3, 2, 1, 0),
             branch = c(F, F, TRUE, F, F, F, F, F, F, F))
)

ttl <- c("Corynebacterium\nAsp -> AMP (purA-purB module)",
         "Psychrobacter / Microbacterium\nFolate biosynthesis (folE / folA)",
         "Psychrobacter\nSulfur relay -> Moco (mec)")

N[, title := ttl[panel]]

E <- rbind(
  data.table(panel = 1, x1 = c(-1.6, -1.6, 0, 0, 0, 0, 0), y1 = c(4, 4, 4, 3, 3, 2, 1),
             x2 = c(-1.6, 0, 0, 1.6, 0, 0, 0), y2 = c(3, 4, 3, 3, 2, 1, 0),
             lab = c("alaT/ald", "TCA", "aspB/GOT", "asnB", "purA [AMG]", "purB [AMG]", " "),
             type = c("Bacterial", "Downstream", "Bacterial", "Bacterial", "Bacteria+Viral", "Bacteria+Viral", "Downstream"),
             flux = c("base", "down", "base", "base", "virus", "virus", "down")),
  data.table(panel = 2, x1 = 0, y1 = 7:1, x2 = 0, y2 = 6:0,
             lab = c("folE [AMG]", "folB", "folK", "folP", "folC", "folA [AMG]", " "),
             type = c("Bacteria+Viral", "Bacterial", "Bacterial", "Bacterial", "Bacterial", "Bacteria+Viral", "Downstream"),
             flux = c("virus", "base", "base", "base", "base", "virus", "down")),
  data.table(panel = 3, x1 = c(0, 0, 0, 0, 0, 2, 2, 2, 2), y1 = c(4, 3, 3, 2, 3, 4, 3, 2, 1),
             x2 = c(0, -1.8, 0, 0, 2, 2, 2, 2, 2), y2 = c(3, 3, 2, 1, 3, 3, 2, 1, 0),
             lab = c("iscS", " ", "carrier-loading", "mec [AMG]", "S relay", "moaA/moaC", "moaD/moaE/moaX", "moeA", " "),
             type = c("Bacterial", "Bacterial", "Bacterial", "Viral AMG", "Bacterial", "Bacterial", "Bacterial", "Bacterial", "Downstream"),
             flux = c("base", "base", "base", "virus", "base", "base", "base", "base", "down"))
)

E[, viral := type %in% c("Viral AMG", "Bacteria+Viral")]
E[, vert := (x1 == x2)]
E[, xs := ifelse(vert, x1, x1 + 0.55 * sign(x2 - x1))]
E[, xe := ifelse(vert, x2, x2 - 0.55 * sign(x2 - x1))]
E[, ys := ifelse(vert, y1 - 0.4 * sign(y1 - y2), y1)]
E[, ye := ifelse(vert, y2 + 0.4 * sign(y1 - y2), y2)]
E[, lx := ifelse(vert, (x1 + x2) / 2 + 0.5, (x1 + x2) / 2)]
E[, ly := ifelse(vert, (y1 + y2) / 2, (y1 + y2) / 2 + 0.3)]
E[, tx := ifelse(vert, (x1 + x2) / 2 - 0.35, (x1 + x2) / 2)]
E[, ty := ifelse(vert, (y1 + y2) / 2, (y1 + y2) / 2 + 0.5)]

cols <- c("Bacterial" = "#2166AC", "Viral AMG" = "#D73027", "Bacteria+Viral" = "#8856A7",
          "Downstream" = "grey60", "Missing" = "grey75")
lts <- c("Bacterial" = "solid", "Viral AMG" = "solid", "Bacteria+Viral" = "solid",
         "Downstream" = "solid", "Missing" = "dashed")

xlims <- list(c(-2.6, 2.6), c(-1.3, 1.9), c(-2.8, 3.4))
ylims <- list(c(-0.8, 4.8), c(-0.8, 7.8), c(-0.8, 4.8))

draw_panel <- function(p) {
  ggplot() +
    geom_segment(data = E[panel == p], aes(x = xs, y = ys, xend = xe, yend = ye, 
                                           color = type, linetype = type, linewidth = flux),
                 arrow = arrow(length = unit(2.6, "mm"))) +
    geom_point(data = E[panel == p & viral == TRUE], aes(x = tx, y = ty), 
               shape = 17, color = "#D73027", size = 2.2) +
    geom_text(data = E[panel == p & lab != " "], aes(x = lx, y = ly, label = lab), 
              size = 3, fontface = "italic", color = "grey20") +
    geom_label(data = N[panel == p], aes(x = x, y = y, label = name, fill = factor(branch)), 
               color = "grey25", size = 2.9) +
    scale_fill_manual(values = c("FALSE" = "grey93", "TRUE" = "#FDB863"), guide = "none") +
    scale_linewidth_manual(values = c(base = 0.7, down = 0.8, virus = 1.6), guide = "none") +
    coord_cartesian(xlim = xlims[[p]], ylim = ylims[[p]]) +
    scale_color_manual(values = cols) + 
    scale_linetype_manual(values = lts) +
    ggtitle(ttl[p]) + 
    theme_void() +
    theme(plot.title = element_text(size = 10.5, face = "bold.italic", hjust = 0.5),
          panel.background = element_rect(fill = "white", color = "grey80"))
}

fig <- draw_panel(1) + draw_panel(2) + draw_panel(3) +
  plot_layout(ncol = 3, widths = c(1, 0.75, 1.2), guides = "collect") +
  plot_annotation(
    title = "Viral AMG-mediated metabolic redirection in early-life colonizers",
    caption = "Thin blue = bacterial baseline flux; thick red/purple = viral-redirected flux; orange box = branch point.\nMetabolite names abbreviated; full KEGG names/reactions see Supplementary Table.",
    theme = theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      plot.caption = element_text(size = 8, hjust = 0.5)
    )
  ) &
  theme(legend.position = "bottom")

ggsave("../output/Fig6/Figure6C_AMG_Schematic.pdf", fig, width = 12, height = 8)
message("✓ Figure 6C schematic saved")
message("✓ All Fig6C analyses completed")
