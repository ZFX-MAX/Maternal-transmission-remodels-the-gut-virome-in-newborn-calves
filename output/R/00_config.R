# ============================================================
# 全局配置：颜色映射、比较列表、通用参数
# ============================================================

# 颜色配置
group_color <- c(
  M21dR = "#EDBA83", M0dR = "#FB9C1D", M21dF = "#9ED5D8", M0dF = "#75CED9",
  C0dME = "#FF9B80", C0dOR = "#5C9E78", C70dF = "#A7C3D9",
  Baoji_C0dOR = "#5C9E78", Baoji_C0dME = "#FF9B80",
  Baoji_M0dF = "#75CED9", Baoji_M21dF = "#9ED5D8",
  Baoji_M0dR = "#FB9C1D", Baoji_M21dR = "#EDBA83",
  Gansu_C0dME = "#FF9B80", Gansu_C70dF = "#A7C3D9",
  Gansu_M0dF = "#75CED9", Gansu_M0dR = "#FB9C1D"
)

# 目标分组
TARGET_GROUPS <- c(
  "Baoji_C0dME", "Baoji_C0dOR", "Baoji_M0dF", "Baoji_M21dF",
  "Baoji_M0dR", "Baoji_M21dR",
  "Gansu_C0dME", "Gansu_C70dF", "Gansu_M0dF", "Gansu_M0dR"
)

# 显著性比较列表（按Area分组）
all_comparisons <- list(
  Baoji = list(
    c("M0dF", "C0dME"), c("M0dR", "C0dME"),
    c("M0dF", "C0dOR"), c("M0dR", "C0dOR")
  ),
  Gansu = list(
    c("C0dME", "C70dF"), c("C70dF", "M0dF"),
    c("M0dF", "M0dR"), c("C70dF", "M0dR"),
    c("C0dME", "M0dF"), c("C0dME", "M0dR")
  )
)
comparisons <- list(c("Baoji_M21dR", "Baoji_M21dF"), c("Baoji_M0dR", "Baoji_M0dF"),
                    c("Baoji_M21dR", "Baoji_M0dR"), c("Baoji_M0dF", "Baoji_M21dF"),
                    c("Baoji_ES", "Baoji_EW"), c("Baoji_CME", "Baoji_M0dF"),c("Baoji_CME", "Baoji_M0dR"),
                    c("Baoji_CME", "Baoji_C0dOR"), c("Baoji_M0dF", "Baoji_C0dOR"),
                    c("Baoji_M0dR", "Baoji_C0dOR"), c("Baoji_ES", "Baoji_M0dR"),
                    c("Gansu_CME", "Gansu_CME70"), c("Gansu_CME70", "Gansu_M0dF"),
                    c("Gansu_M0dF", "Gansu_M0dR"), c("Gansu_CME70", "Gansu_M0dR"),c("Gansu_CME", "Gansu_M0dR") )
# 核心属列表（用于韦恩图和下游分析）
target_genus <- c(
  "Agathobacter", "Aerococcus", "Anaerostipes", "Arthrobacter",
  "Bacteroides", "Blautia", "Brachybacterium", "Corynebacterium",
  "Dorea", "Enterocloster", "Escherichia", "Faecalibacterium",
  "Glutamicibacter", "Holdemanella", "Jeotgalicoccus",
  "Mediterraneibacter", "Megasphaera", "Microbacterium",
  "Moraxella", "Parabacteroides", "Phocaeicola", "Pseudoalteromonas",
  "Pseudomonas", "Psychrobacter", "Rhodococcus", "Staphylococcus",
  "Stenotrophomonas", "Streptomyces", "Subdoligranulum"
)

# 母系传播路径定义
MOTHER_GROUPS <- c(
  "C0dOR_M0dR", "C0dOR_M21dR", "C0dOR_M21dF", "C0dOR_M0dF",
  "C0dME_M0dR", "C0dME_M0dF", "C0dME_M21dR", "C0dME_M21dF",
  "C70dF_M0dR", "C70dF_M0dF"
)
theme_pub <- function(base_size = 12, base_family = "") {
  theme_classic(base_size = base_size, base_family = base_family) %+replace%
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = rel(1.1)),
      axis.title = element_text(size = rel(1.0)),
      axis.text  = element_text(size = rel(0.9)),
      legend.title = element_text(face = "bold", size = rel(0.95)),
      legend.text  = element_text(size = rel(0.85)),
      strip.text = element_text(face = "bold", size = rel(0.9)),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}
message("✓ Global configuration loaded")
