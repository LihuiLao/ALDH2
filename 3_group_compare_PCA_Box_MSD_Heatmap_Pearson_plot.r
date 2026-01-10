# 1.  (Load Libraries)
# 
# install.packages(c("ggplot2", "readxl", "psych", "pheatmap", "limma", "reshape2", "writexl", "dplyr"))
library(ggplot2)
library(readxl)
library(psych)
library(pheatmap)
library(limma)
library(reshape2) #  melt 
library(writexl)  # ：Excel
library(dplyr)    # ：

# 2.  (Define Column Names and Parameters)
# 
mut_cols <- c("Mut_1", "Mut_2", "Mut_3")
wt_cols <- c("WT_1", "WT_2", "WT_3")
ctl_cols <- c("Con_1", "Con_2", "Con_3")

# 
intensity_cols <- c(mut_cols, wt_cols, ctl_cols)

# ANOVA，

# 3.  (Data Loading and Pre-processing)
#  'lfq_analysis' sheet 9
df <- read_xlsx('analysis.xlsx', sheet = 'lfq_analysis')

# QC
# ，df2
df2 <- df[, intensity_cols]

# 4.  (Data Visualization)

# 4.1  (Boxplot of Sample Intensities)
df2m <- reshape2::melt(df2, measure.vars = intensity_cols)
color_map <- setNames(
  c(rep("red", length(mut_cols)),
    rep("green", length(wt_cols)),
    rep("blue", length(ctl_cols))),
  intensity_cols
)
ggplot(df2m, aes(x = variable, y = value, fill = variable)) +
  geom_boxplot() +
  scale_fill_manual(values = color_map) + 
  theme_bw() +
  theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "", x = 'Sample', y = 'Log2 LFQ Intensity')

# 4.2  (Pairs Panels)
pairs.panels(df2)

# 4.3  (MDS Plot)
mds_colors <- ifelse(colnames(df2) %in% mut_cols, "red",
                      ifelse(colnames(df2) %in% wt_cols, "green", "blue"))
plotMDS(df2, col = mds_colors, labels = colnames(df2))
legend("topleft", legend = c("Mutant", "Wild-Type", "Control"), fill = c("red", "green", "blue"))

# 4.4  (PCA Plot)
pca_data <- t(df2)
pca_result <- prcomp(pca_data, scale. = TRUE)
pca_plot_data <- as.data.frame(pca_result$x)
pca_plot_data$group <- factor(
  ifelse(rownames(pca_plot_data) %in% mut_cols, "Mutant",
         ifelse(rownames(pca_plot_data) %in% wt_cols, "Wild-Type", "Control")),
  levels = c("Mutant", "Wild-Type", "Control")
)
percent_var <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100
ggplot(pca_plot_data, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 5, alpha = 0.8) +
  geom_text(aes(label = rownames(pca_plot_data)), vjust = -1.5, color = "black", size = 3) +
  theme_bw(base_size = 14) +
  scale_color_manual(values = c("Mutant" = "red", "Wild-Type" = "green", "Control" = "blue")) +
  labs(
    title = " (PCA)",
    x = paste0("PC1 (", round(percent_var[1], 1), "% variance)"),
    y = paste0("PC2 (", round(percent_var[2], 1), "% variance)"),
    color = "Group"
  ) +
  theme(legend.position = "bottom")

#  ()
#  intensity_cols 
group_vector <- rep(c("Mutant", "Wild-Type", "Control"), c(length(mut_cols), length(wt_cols), length(ctl_cols)))
names(group_vector) <- c(mut_cols, wt_cols, ctl_cols)
group_vector <- group_vector[intensity_cols] #  intensity_cols 

annotation_col <- data.frame(Group = factor(group_vector, levels = c("Mutant", "Wild-Type", "Control")))
rownames(annotation_col) <- intensity_cols
annotation_colors <- list(Group = c(Mutant = "red", `Wild-Type` = "green", Control = "blue"))


# 4.5  ()
#  pheatmap 
heatmap_all_clustered <- pheatmap(df2,
         scale = 'row',
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         cluster_cols = TRUE,               # <-- 
         cluster_rows = TRUE,
         show_rownames = FALSE,
         main = "/ ()",
         silent = FALSE) #  silent=TRUE 

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
         cluster_cols = FALSE,              # <-- 
         cluster_rows = TRUE,
         show_rownames = FALSE,
         main = "/ ()")


# 5.  (Specific Gene Set Analysis)
dg <- read_xlsx('analysis.xlsx', sheet = 'fatty_acid_related')
if (all(intensity_cols %in% colnames(dg))) {
  dg_intensity <- dg[, intensity_cols]

  # 5.1  ()
  #  pheatmap 
  heatmap_fa_clustered <- pheatmap(dg_intensity,
           scale = 'row',
           annotation_col = annotation_col,
           annotation_colors = annotation_colors,
           cluster_cols = TRUE,               # <-- 
           cluster_rows = TRUE,
           labels_row = dg$`Leading_gene`,
           main = " ()",
           silent = FALSE) #  silent=TRUE 

  # 
  row_order_fa <- heatmap_fa_clustered$tree_row$order
  #  'dg' 
  dg_heatmap_ordered <- dg[row_order_fa, ]
  # Excel
  write_xlsx(dg_heatmap_ordered, "heatmap_fatty_acid_ordered.xlsx")
  # ，
  print("“” heatmap_fatty_acid_ordered.xlsx ")

  # 5.2  ()
  pheatmap(dg_intensity,
           scale = 'row',
           annotation_col = annotation_col,
           annotation_colors = annotation_colors,
           cluster_cols = FALSE,              # <-- 
           cluster_rows = TRUE,
           labels_row = dg$`Leading_gene`,
           main = " ()")

} else {
  print(": 'fatty_acid_related' ，")
}
