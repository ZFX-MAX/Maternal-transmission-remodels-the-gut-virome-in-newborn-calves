# ============================================================
# Figure 6A-6B: LEfSe Differential Abundance Analysis (Pathway Only)
# Description: Identify differentially enriched viral/bacterial 
#              KEGG pathways across developmental stages
# ============================================================
suppressPackageStartupMessages({
  library(microeco)
  library(data.table)
})

source("../R/00_config.R")

dir.create("../output/Fig6", recursive = TRUE, showWarnings = FALSE)

# ==================== 1. Configuration ====================
GROUP_CONFIG <- list(
  Baoji = c("Baoji_C0dOR", "Baoji_C0dME", "Baoji_M0dF", "Baoji_M0dR"),
  Gansu = c("Gansu_C0dME", "Gansu_C70dF", "Gansu_M0dF")
)

INPUT_FILES <- list(
  virus_pathway = "../data/kegg/AMG_Pathway_abundance.tsv",
  bact_pathway = "../data/kegg/kegg_Pathway_abundance.tsv",
  meta = "../data/meta.tsv"
)

# ==================== 2. LEfSe Analysis Function ====================
# ============================================================
# Figure 6A-6B: LEfSe Differential Abundance Analysis (Pathway Only)
# Description: Identify differentially enriched viral/bacterial 
#              KEGG pathways across developmental stages
# ============================================================
suppressPackageStartupMessages({
  library(microeco)
  library(data.table)
})

source("../R/00_config.R")

dir.create("../output/Fig6", recursive = TRUE, showWarnings = FALSE)

# ==================== 1. Configuration ====================
GROUP_CONFIG <- list(
  Baoji = c("Baoji_C0dOR", "Baoji_C0dME", "Baoji_M0dF", "Baoji_M0dR"),
  Gansu = c("Gansu_C0dME", "Gansu_C70dF", "Gansu_M0dF")
)

INPUT_FILES <- list(
  virus_pathway = "../data/kegg/AMG_Pathway_abundance.tsv",
  bact_pathway = "../data/kegg/kegg_Pathway_abundance.tsv",
  meta = "../data/meta.tsv"
)

# ==================== 2. LEfSe Analysis Function ====================
run_lefse_pathway <- function(abundance_file, meta_file, group_vector, 
                              organism_type, output_dir) {
  message(sprintf("Running LEfSe for %s pathway...", organism_type))
  
  # 读取丰度矩阵
  abundance <- fread(abundance_file)
  feature_col <- colnames(abundance)[1]
  
  pathway_names <- as.character(abundance[[feature_col]])
  abundance[, (feature_col) := NULL]
  abundance <- as.data.frame(abundance)
  rownames(abundance) <- pathway_names
  
  pathway_table <- data.frame(
    pathway = pathway_names,
    row.names = pathway_names,
    stringsAsFactors = FALSE
  )
  
  meta <- fread(meta_file)
  meta <- meta[Group_Area %in% group_vector, ]
  
  group_data <- data.frame(
    Sample_ID = meta$Sample_ID,
    group = meta$Group_Area,
    row.names = meta$Sample_ID,
    stringsAsFactors = FALSE
  )
  
  common_samples <- intersect(rownames(group_data), colnames(abundance))
  group_data <- group_data[common_samples, , drop = FALSE]
  abundance <- abundance[, common_samples, drop = FALSE]
  
  if (nrow(group_data) == 0) {
    message("  Warning: No matching samples found!")
    return(data.frame())
  }
  
  mt <- microtable$new(
    sample_table = group_data,
    otu_table = abundance,
    tax_table = pathway_table
  )
  
  set.seed(1234)
  lefse_res <- trans_diff$new(
    dataset = mt,
    method = "lefse",
    group = "group",
    taxa_level = "pathway"
  )
  
  diff_res <- as.data.frame(lefse_res$res_diff)
  diff_res_sig <- diff_res[diff_res$LDA > 3 & diff_res$P.adj < 0.01, ]
  if ("Taxa" %in% colnames(diff_res_sig)) {
    n_before <- nrow(diff_res_sig)
    diff_res_sig <- diff_res_sig[!grepl("^map", diff_res_sig$Taxa), ]
    message(sprintf("  Filtered out %d 'map'-prefixed duplicate pathways", 
                    n_before - nrow(diff_res_sig)))
  }
  
  cohort_name <- if (any(grepl("Baoji", group_vector))) "Baoji" else "Gansu"
  output_file <- file.path(
    output_dir,
    sprintf("%s_%s_lefse_pathway.csv", cohort_name, organism_type)
  )
  fwrite(diff_res_sig, output_file)
  
  message(sprintf("  Found %d significant pathways, saved to %s", 
                  nrow(diff_res_sig), output_file))
  
  return(diff_res_sig)
}

# ==================== 3. Run LEfSe Analyses (Pathway Only) ====================
message("Starting LEfSe pathway analyses...")

for (org in c("virus", "bact")) {
  run_lefse_pathway(
    abundance_file = INPUT_FILES[[paste0(org, "_pathway")]],
    meta_file = INPUT_FILES$meta,
    group_vector = GROUP_CONFIG$Baoji,
    organism_type = org,
    output_dir = "../output/Fig6"
  )
}

for (org in c("virus", "bact")) {
  run_lefse_pathway(
    abundance_file = INPUT_FILES[[paste0(org, "_pathway")]],
    meta_file = INPUT_FILES$meta,
    group_vector = GROUP_CONFIG$Gansu,
    organism_type = org,
    output_dir = "../output/Fig6"
  )
}

message("✓ LEfSe pathway analyses completed")
