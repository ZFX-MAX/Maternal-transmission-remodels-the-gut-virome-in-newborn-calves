# ============================================================
# Figure 4A-4B: 共享事件鉴定与韦恩图分析
# ============================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggvenn)
  library(patchwork)
})

source("../R/00_config.R")

# 创建输出目录
dir.create("../output/Fig4", recursive = TRUE, showWarnings = FALSE)

# ==================== 1. 数据加载与处理 ====================
PATHS <- list(
  vir1   = "../data/cluster_viruse.tsv",
  vir2   = "..//data/cluster_no_0dKG.tsv",
  bact   = "../data/cluster_bact.tsv",
  anno_v = "../data/blast_viruse4.tsv",
  anno_b = "../data/blast_bact.tsv",
  meta   = "../data/meta.tsv"
)

process_cluster <- function(cluster_file, anno_file, meta_file, type) {
  dt <- fread(cluster_file, header = FALSE, sep = "\t",
              col.names = c("Cluster", "Contigs"))
  valid <- dt[, .N, by = Cluster][N >= 3, Cluster]
  dt <- dt[Cluster %in% valid & Contigs != Cluster]
  
  dt[, Cluster_Prefix := sub("\\.NODE_.*", "", Cluster)]
  dt[, Contig_Prefix  := sub("\\.NODE_.*", "", Contigs)]
  
  meta <- fread(meta_file, header = TRUE)
  meta[, Sample_ID := as.character(Sample_ID)]
  meta_core <- unique(meta[, .(Sample_ID, Individuality, Paired_number,
                               Age, Group, Area)])
  
  mc <- copy(meta_core)
  setnames(mc, c("Sample_ID","Individuality","Paired_number","Age","Group","Area"),
           c("Cluster_Prefix","Individuality_cluster","Paired_number_cluster",
             "Age_cluster","Group_cluster","Area_cluster"))
  dt <- mc[dt, on = "Cluster_Prefix"]
  
  mg <- copy(meta_core)
  setnames(mg, c("Sample_ID","Individuality","Paired_number","Age","Group","Area"),
           c("Contig_Prefix","Individuality_Contigs","Paired_number_Contigs",
             "Age_Contigs","Group_Contigs","Area_Contigs"))
  dt <- mg[dt, on = "Contig_Prefix"]
  
  anno <- unique(fread(anno_file), by = "Contig_name")
  if (type == "Viruse") {
    anno_tax <- anno[, .(Contig_name, Annotation = host_annotation2,
                         Predicted_label = predicted_label)]
  } else {
    anno_tax <- anno[, .(Contig_name, Annotation = annotation_pre,
                         Predicted_label = NA)]
  }
  
  ac <- copy(anno_tax)
  setnames(ac, c("Contig_name","Annotation","Predicted_label"),
           c("Cluster","Annotation_cluster","Predicted_label_cluster"))
  dt <- ac[dt, on = "Cluster"]
  dt[, `:=`(Annotation_Contigs = Annotation_cluster,
            Predicted_label_Contigs = Predicted_label_cluster)]
  dt[, `:=`(Sample_cluster = Cluster_Prefix, Sample_Contigs = Contig_Prefix,
            Cluster_Prefix = NULL, Contig_Prefix = NULL)]
  
  std_cols <- c("Contigs","Cluster","Sample_Contigs","Sample_cluster",
                "Annotation_Contigs","Annotation_cluster",
                "Predicted_label_Contigs","Predicted_label_cluster",
                "Individuality_Contigs","Individuality_cluster",
                "Paired_number_Contigs","Paired_number_cluster",
                "Age_Contigs","Age_cluster",
                "Group_Contigs","Group_cluster",
                "Area_Contigs","Area_cluster")
  for (col in std_cols) if (!col %in% names(dt)) dt[, (col) := NA_character_]
  dt <- dt[, ..std_cols]
  dt[, Type := type]
  return(dt)
}

# 构建并过滤数据
vir1 <- process_cluster(PATHS$vir1, PATHS$anno_v, PATHS$meta, "Viruse")
vir2 <- process_cluster(PATHS$vir2, PATHS$anno_v, PATHS$meta, "Viruse")
bact <- process_cluster(PATHS$bact, PATHS$anno_b, PATHS$meta, "Bacteria")

cluster_viruse <- rbindlist(list(vir1, vir2), use.names = TRUE)[!duplicated(Contigs)]
cluster_merge  <- rbindlist(list(bact, cluster_viruse), use.names = TRUE)

target_grps <- c("C0dOR","M0dR","C0dME","C70dF","M0dF","M21dR","M21dF")
cross_pat   <- "C0dME|C0dOR|C70dF"

dt <- cluster_merge[
  Group_cluster %in% target_grps &
    (grepl(cross_pat, Group_Contigs) | grepl(cross_pat, Group_cluster)) &
    Area_cluster == Area_Contigs &
    Group_cluster != Group_Contigs &
    !is.na(Annotation_cluster) & Annotation_cluster != "" &
    !is.na(Group_Contigs) & Group_Contigs != "" &
    Paired_number_Contigs != "Baoji_NA"
]

min_freq <- 10
freq_table <- dt[, .N, by = Annotation_cluster]
dt <- dt[Annotation_cluster %in% freq_table[N >= min_freq, Annotation_cluster]]

dt[, Groups := paste(pmin(Group_cluster, Group_Contigs),
                     pmax(Group_cluster, Group_Contigs), sep = "_")]
dt[, Paired_Type := fifelse(Paired_number_cluster == Paired_number_Contigs,
                            "Paired","Unpaired")]

# 保存中间对象供下游脚本使用
dir.create("../output/intermediate", recursive = TRUE, showWarnings = FALSE)
saveRDS(dt, "../output/intermediate/shared_events_dt.rds")
fwrite(dt, "../output/intermediate/shared_events_dt.tsv",sep="\t")
# ==================== 2. Figure 4A: 柱状图 ====================
plot_dt <- dt[!is.na(Groups) & Groups != "",
              .(Count = .N), by = .(Region = Area_cluster, Type, Groups, Paired_Type)]
plot_dt[, Groups := factor(Groups, levels = sort(unique(Groups)))]
plot_dt[, Paired_Type := factor(Paired_Type, levels = c("Paired","Unpaired"))]

p_4a <- ggplot(plot_dt, aes(x = Groups, y = log10(Count + 1), fill = Paired_Type)) +
  geom_col(position = position_dodge2(preserve = "single"), width = 0.7) +
  facet_grid(Type ~ Region, scales = "free", drop = TRUE, space = "free") +
  scale_y_continuous(
    breaks = c(0, 1, 2, 3, 4),
    labels = c("0", "10", "100", "1000", "10000")
  ) +
  scale_fill_manual(values = c("Paired" = "#7A2323", "Unpaired" = "#BEBEBE")) +
  labs(x = "Groups",
       y = expression("Number of Contigs (" * log[10] * "(x+1))"),
       fill = "Pairing Status") +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold", size = 11),
        strip.background = element_rect(fill = "#F5F5F5", color = NA),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "gray", fill = NA, linewidth = 0.5),
        legend.position = "right")




ggsave("../output/Fig4/Figure4A_Shared_Events.pdf", p_4a, width = 8, height = 5)
message("✓ Figure 4A saved")


plot_dt <- unique(dt[!is.na(Groups) & Groups != "",
                     .(Region=Area_cluster, Type, Groups,
                       Sample1=Sample_cluster, Sample2=Sample_Contigs,
                       Pair1=Paired_number_cluster, Pair2=Paired_number_Contigs,
                       Cluster)])[
                         , .(Shared=uniqueN(Cluster)),
                         by=.(Region,Type,Groups,Sample1,Sample2,Pair1,Pair2)
                       ][
                         , Paired_Type:=ifelse(Pair1==Pair2,"Paired","Unpaired")
                       ]
p_4a <- ggplot(plot_dt,aes(Groups,Shared+1,fill=Paired_Type)) +
  geom_boxplot(position=position_dodge(.75),width=.6,
               outlier.shape=NA,alpha=.65) +
  geom_jitter(aes(color=Paired_Type),
              position=position_jitterdodge(.12,dodge.width=.75),
              size=1.3,alpha=.35) +
  scale_y_log10() +
  facet_grid(Type~Region,scales="free_x",space="free") +
  scale_fill_manual(values=c(Paired="#7A2323",Unpaired="#BEBEBE")) +
  scale_color_manual(values=c(Paired="#7A2323",Unpaired="#888888")) +
  labs(x="Sample groups",y="Shared viral clusters (log scale)",
       fill="Pairing status",color="Pairing status") +
  theme_minimal(base_size=12) +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        panel.grid.major.x=element_blank(),
        panel.grid.minor=element_blank())

# ==================== 3. Figure 4B: 韦恩图 ====================
meta_shared <- rbindlist(list(
  unique(dt[, .(ID = Contigs, Type, Sample = Sample_Contigs,
                Annotation = Annotation_Contigs,
                Predicted_Label = Predicted_label_Contigs,
                Individuality = Individuality_Contigs,
                Paired_Number = Paired_number_Contigs,
                Age = Age_Contigs, Group = Group_Contigs, Area = Area_Contigs)]),
  unique(dt[, .(ID = Cluster, Type, Sample = Sample_cluster,
                Annotation = Annotation_cluster,
                Predicted_Label = Predicted_label_cluster,
                Individuality = Individuality_cluster,
                Paired_Number = Paired_number_cluster,
                Age = Age_cluster, Group = Group_cluster, Area = Area_cluster)])
), use.names = TRUE)
setorder(meta_shared, Type, ID)

cc_map <- unique(dt[, .(Cluster_ID = Cluster, Contig_ID = Contigs)])
all_pairs <- cc_map[, .(Cluster_ID, Event_Type = "Contig_vs_Cluster",
                        Left_Item = Contig_ID, Right_Item = Cluster_ID)]

setkey(meta_shared, ID)
fc <- copy(all_pairs)
fc <- meta_shared[fc, on = c(ID = "Left_Item")]
setnames(fc, setdiff(names(meta_shared),"ID"), paste0("L_", setdiff(names(meta_shared),"ID")))
setnames(fc, "ID", "Left_Item")
fc <- meta_shared[fc, on = c(ID = "Right_Item")]
setnames(fc, setdiff(names(meta_shared),"ID"), paste0("R_", setdiff(names(meta_shared),"ID")))
setnames(fc, "ID", "Right_Item")

fc2 <- fc[L_Group != R_Group]
fc2[, Paired_Type := fifelse(L_Paired_Number == R_Paired_Number, "Paired","Unpaired")]
fc2[, Groups := {
  g1 <- fifelse(is.na(L_Group)|L_Group=="", NA_character_, L_Group)
  g2 <- fifelse(is.na(R_Group)|R_Group=="", NA_character_, R_Group)
  fifelse(is.na(g1)&is.na(g2), NA_character_,
          fifelse(is.na(g1), g2, fifelse(is.na(g2), g1,
                                         paste(pmin(g1,g2), pmax(g1,g2), sep="_"))))
}]

cond_list <- list(
  "Calf oral to Meconium (Paired)"   = fc2[Groups=="C0dME_C0dOR" & Paired_Type=="Paired"],
  "Calf oral to Meconium (Unpaired)" = fc2[Groups=="C0dME_C0dOR" & Paired_Type=="Unpaired"],
  "Mother to Offspring (Paired)"     = fc2[Groups %in% MOTHER_GROUPS & Paired_Type=="Paired"],
  "Mother to Offspring (Unpaired)"   = fc2[Groups %in% MOTHER_GROUPS & Paired_Type=="Unpaired"]
)
target_genus <- fread("../output/Fig2/LEFse_overlap.tsv",header = F,sep=" ")
target_genus <- rbind(target_genus,fread("../output/Fig3/LEFse_overlap.tsv",header = F,sep=" ") )
target_genus <- target_genus[target_genus$V2%in%c("Gansu_C0dME","Gansu_C70dF","Baoji_C0dOR","Baoji_C0dME"),]
target_genus <- unique(unlist(strsplit(target_genus$V6, split = ";")))
process_venn <- function(name, dt_sub) {
  anno <- dt_sub$L_Annotation
  type <- dt_sub$L_Type
  T_set <- target_genus
  B_set <- unique(anno[!is.na(type) & type == "Bacteria"])
  V_set <- unique(anno[!is.na(type) & type == "Viruse"])
  
  sets <- list("Offspring Enriched" = T_set,
               "Shared Bacteria" = B_set,
               "Shared Viruses" = V_set)
  
  regions <- list(
    "A" = setdiff(T_set, union(B_set, V_set)),
    "B" = setdiff(B_set, union(T_set, V_set)),
    "C" = setdiff(V_set, union(T_set, B_set)),
    "D" = setdiff(intersect(T_set, B_set), V_set),
    "E" = setdiff(intersect(T_set, V_set), B_set),
    "F" = setdiff(intersect(B_set, V_set), T_set),
    "G" = Reduce(intersect, list(T_set, B_set, V_set)))
  
  df_out <- data.frame(Condition = name,
                       Venn_Label = rep(names(regions), lengths(regions)),
                       Genera = unlist(regions), stringsAsFactors = FALSE)
  
  p <- ggvenn(sets, fill_color = c("#999999","#56B4E9","#E69F00"),
              fill_alpha = 0.5, stroke_color = "black",
              set_name_size = 3, text_size = 3,
              show_percentage = FALSE, stroke_size = 0.5, padding = 0.1) +
    ggtitle(name) + theme_void() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10))
  
  list(plot = p, table = df_out)
}

res_list <- lapply(names(cond_list), function(nm) process_venn(nm, cond_list[[nm]]))
final_plot <- (res_list[[1]]$plot + res_list[[2]]$plot) /
  (res_list[[3]]$plot + res_list[[4]]$plot)

final_table <- do.call(rbind, lapply(res_list, function(x) x$table))
final_table$Meaning <- factor(final_table$Venn_Label,
                              levels = c("A","B","C","D","E","F","G"),
                              labels = c("Concordant bacterial and phage-host enrichment only",
                                         "Bacterial sequence sharing only",
                                         "Predicted host genera of shared phages only",
                                         "Concordant bacterial and phage-host enrichment + bacterial sequence sharing",
                                         "Concordant bacterial and phage-host enrichment + predicted host genera of shared phages",
                                         "Bacterial sequence sharing + predicted host genera of shared phages",
                                         "Bacterial enrichment + phage-host enrichment + bacterial sequence sharing + predicted host genera of shared phages"))


# 保存final_table供下游使用
saveRDS(final_table, "../output/intermediate/final_table.rds")
fwrite(final_table, "../output/intermediate/final_table.txt",sep="\t")
ggsave("../output/Fig4/Figure4B_Venn.pdf", final_plot, width = 4.5, height = 4.5)
message("✓ Figure 4B saved")
message("✓ Shared events analysis completed")

