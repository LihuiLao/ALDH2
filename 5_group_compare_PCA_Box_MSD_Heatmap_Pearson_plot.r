# 1.  (Load Libraries)
# 
# install.packages(c("ggplot2", "readxl", "psych", "pheatmap", "limma", "reshape2", "writexl", "dplyr", "RColorBrewer"))
library(ggplot2)
library(readxl)
library(psych)
library(pheatmap)
library(limma)
library(reshape2) #  melt 
library(writexl)  # Excel
library(dplyr)    # 
library(RColorBrewer) # 

# 2.  (Define Column Names and Parameters)
# ---  ---
# 
# 
mut_cols <- c("LLH1_1", "LLH1_2", "LLH1_3")
wt_cols  <- c("LLH10_1", "LLH10_2", "LLH10_3")
ctl_cols <- c("LLH24_1", "LLH24_2", "LLH24_3")
ng_cols  <- c("LLH4_1", "LLH4_2", "LLH4_3")
group5_cols <- c("Group5_Sample1", "Group5_Sample2", "Group5_Sample3") 

# ，
group_list <- list(
  "Mutant" = mut_cols,
  "Wild-Type" = wt_cols,
  "Control" = ctl_cols,
  "NewGroup" = ng_cols,
  "Group5" = group5_cols  # <-- 
)

#  (，)
intensity_cols <- unlist(group_list)

#  (，)
group_names <- names(group_list)
# RColorBrewer，
# colors <- c("red", "green", "blue", "purple", "orange")
colors <- brewer.pal(length(group_names), "Set1")
names(colors) <- group_names


# 3.  (Data Loading and Pre-processing)
# ---  ---
#  'analysis.xlsx'  'lfq_analysis' sheet 
df <- read_xlsx('analysis.xlsx', sheet = 'lfq_analysis')

# QC
# ，df2
#  intensity_cols  df 
if (!all(intensity_cols %in% colnames(df))) {
  stop(": Excel，2 'group_list' Excel")
}
df2 <- df[, intensity_cols]


# 4.  (Data Visualization)
# ，

# 4.1  (Boxplot of Sample Intensities)
df2m <- reshape2::melt(df2, measure.vars = intensity_cols)
# 
color_map_boxplot <- unlist(sapply(group_names, function(g) rep(colors[g], length(group_list[[g]]))))
names(color_map_boxplot) <- intensity_cols

ggplot(df2m, aes(x = variable, y = value, fill = variable)) +
  geom_boxplot() +
  scale_fill_manual(values = color_map_boxplot) +
  theme_bw() +
  theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "", x = 'Sample', y = 'Log2 LFQ Intensity')

# 4.2  (Pairs Panels)
# ：，，
pairs.panels(df2)

# 4.3  (MDS Plot)
# MDS
group_membership <- sapply(colnames(df2), function(col) {
  names(which(sapply(group_list, function(g) col %in% g)))
})
mds_colors <- colors[group_membership]

plotMDS(df2, col = mds_colors, labels = colnames(df2))
legend("topleft", legend = group_names, fill = colors[group_names])

# 4.4  (PCA Plot)
pca_data <- t(df2)
pca_result <- prcomp(pca_data, scale. = TRUE)
pca_plot_data <- as.data.frame(pca_result$x)

# 
pca_plot_data$group <- factor(group_membership, levels = group_names)

percent_var <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100
ggplot(pca_plot_data, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 5, alpha = 0.8) +
  geom_text(aes(label = rownames(pca_plot_data)), vjust = -1.5, color = "black", size = 3) +
  theme_bw(base_size = 14) +
  scale_color_manual(values = colors) +
  labs(
    title = " (PCA)",
    x = paste0("PC1 (", round(percent_var[1], 1), "% variance)"),
    y = paste0("PC2 (", round(percent_var[2], 1), "% variance)"),
    color = "Group"
  ) +
  theme(legend.position = "bottom")

#  ()
group_vector <- factor(group_membership, levels = group_names)
annotation_col <- data.frame(Group = group_vector)
rownames(annotation_col) <- intensity_cols
annotation_colors <- list(Group = colors)


# 4.5  ()
#  pheatmap 
heatmap_all_clustered <- pheatmap(df2,
         scale = 'row',
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         cluster_cols = TRUE,          # <-- 
         cluster_rows = TRUE,
         show_rownames = FALSE,
         main = "/ ()",
         silent = TRUE) #  silent=TRUE 

# 
row_order_all <- heatmap_all_clustered$tree_row$order
#  'df'
df_heatmap_ordered <- df[row_order_all, ]
# Excel
write_xlsx(df_heatmap_ordered, "heatmap_all_proteins_ordered.xlsx")
# ，
print("“” heatmap_all_proteins_ordered.xlsx ")


# 4.6  ()
pheatmap(df2,
         scale = 'row',
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         cluster_cols = FALSE,         # <-- 
         cluster_rows = TRUE,
         show_rownames = FALSE,
         main = "/ ()")


# 5.  (Specific Gene Set Analysis)
#  'fatty_acid_related' 
# ---  ---
# ， 'fatty_acid_related' sheet 
dg <- read_xlsx('analysis.xlsx', sheet = 'fatty_acid_related')

# 
if (all(intensity_cols %in% colnames(dg))) {
  dg_intensity <- dg[, intensity_cols]

  # 5.1  ()
  #  pheatmap 
  heatmap_specific_clustered <- pheatmap(dg_intensity,
           scale = 'row',
           annotation_col = annotation_col,
           annotation_colors = annotation_colors,
           cluster_cols = TRUE,          # <-- 
           cluster_rows = TRUE,
           labels_row = dg$`Leading_gene`, # dg
           main = " ()",
           silent = TRUE) #  silent=TRUE 

  # 
  row_order_specific <- heatmap_specific_clustered$tree_row$order
  #  'dg' 
  dg_heatmap_ordered <- dg[row_order_specific, ]
  # Excel
  write_xlsx(dg_heatmap_ordered, "heatmap_specific_genes_ordered.xlsx")
  # ，
  print("“” heatmap_specific_genes_ordered.xlsx ")

  # 5.2  ()
  pheatmap(dg_intensity,
           scale = 'row',
           annotation_col = annotation_col,
           annotation_colors = annotation_colors,
           cluster_cols = FALSE,         # <-- 
           cluster_rows = TRUE,
           labels_row = dg$`Leading_gene`,
           main = " ()")

} else {
  print(": 'fatty_acid_related' ，")
}
