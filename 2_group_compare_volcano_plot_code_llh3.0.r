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
library(dplyr)    # ：，

# 2.  (Define Column Names and Parameters)
# 
mut_cols <- c("Mut_1", "Mut_2", "Mut_3")
wt_cols <- c("WT_1", "WT_2", "WT_3")
# 
intensity_cols <- c(mut_cols, wt_cols)
# 
fc_col <- "FC"
# p
pval_col <- "p-value"

# 
s0 <- 0.1      # Fudge factor，
f <- 5         #  (degrees of freedom)  = -+1
confidence <- 0.95 # 
ta <- qt(confidence, f) # t

# 3.  (Data Loading and Pre-processing)
#  'analysis.xlsx' R，
df <- read_xlsx('analysis.xlsx', sheet = 'lfq_analysis')

# 
data_calsmoothcurve <- function(x, ta, s0, f) {
  y <- ta * (1 + (s0 / ((abs(x) / ta) - s0)))
  y <- -log10(2 * (1 - pt(y, df = f)))
  return(y)
}

# -log10(p-value)
df$lgpvalue <- -log10(df[[pval_col]])

# y
df$fudge <- apply(df, 1, function(x) {
  data_calsmoothcurve(as.numeric(x[fc_col]), ta, s0, f)
})

# 
# 'sig'，'0'
df$sig <- '0'

# （mask）
# 1: p-value (lgpvalue > fudge)
# 2: FC (FC > ta * s0  FC < -ta * s0)

# / ('1')
up_regulated_mask <- (df$lgpvalue > df$fudge) & (df[[fc_col]] > ta * s0)
df$sig[up_regulated_mask] <- '1'

# / ('-1')
down_regulated_mask <- (df$lgpvalue > df$fudge) & (df[[fc_col]] < -ta * s0)
df$sig[down_regulated_mask] <- '-1'

# ()sig，
df$sig <- factor(df$sig, levels = c('-1', '0', '1'))


# 4. Excel (Export Volcano Plot Data to Excel)
#  'Protein IDs'  'Gene names' ，
# ， select() 
df_to_export <- df %>%
  select(any_of(c("Protein IDs", "Gene names")), !!fc_col, !!pval_col, lgpvalue, sig, fudge, everything())

# （lgpvalue, fudge, sig）Excel
# R
write_xlsx(df_to_export, "volcano_analysis_output.xlsx")

# ，
print(" volcano_analysis_output.xlsx ")


# 5.  (Data Visualization)

# 5.1  (Volcano Plot)
ggplot(df, aes_string(x = fc_col, y = "lgpvalue", color = "sig")) +
  geom_point(size = 2, shape = 1, stroke = 1) +
  ylim(0, 8) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("-1" = "green", "0" = "gray", "1" = "red"),
                     name = "Significance") +
  geom_function(fun = function(x) { data_calsmoothcurve(x, ta, s0, f) },
                xlim = c(0.02, 8), linetype = 'dashed', n = 10000,
                linewidth = 1, alpha = 0.3, color = 'darkred') +
  geom_function(fun = function(x) { data_calsmoothcurve(x, ta, s0, f) },
                xlim = c(-8, -0.02),
                linewidth = 1, alpha = 0.3, color = 'darkgreen',
                n = 10000, linetype = 'dashed') +
  labs(x = 'Log2 Fold Change (Mut - WT)', y = '-Log10(P-value)')


# 5.2  (Boxplot of Sample Intensities)
df2 <- df[, intensity_cols]
df2m <- reshape2::melt(df2, measure.vars = intensity_cols)
color_map <- setNames(
  c(rep("red", length(mut_cols)), rep("green", length(wt_cols))),
  intensity_cols
)
ggplot(df2m, aes(x = variable, y = value, fill = variable)) +
  geom_boxplot() +
  scale_fill_manual(values = color_map) +
  theme_bw() +
  theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = 'Sample', y = 'Log2 LFQ Intensity')


# 5.3  (Pairs Panels for Sample Correlation)
pairs.panels(df2)


# 5.4  (MDS Plot)
mds_colors <- ifelse(colnames(df2) %in% mut_cols, "red", "green")
plotMDS(df2, col = mds_colors, labels = colnames(df2))
legend("topleft", legend = c("Mutant", "Wild-Type"), fill = c("red", "green"))

# 5.5  (PCA Plot)
pca_data <- t(df2)
pca_result <- prcomp(pca_data, scale. = TRUE)
pca_plot_data <- as.data.frame(pca_result$x)
pca_plot_data$group <- ifelse(rownames(pca_plot_data) %in% mut_cols, "Mutant", "Wild-Type")
percent_var <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100
ggplot(pca_plot_data, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 5, alpha = 0.8) +
  geom_text(aes(label = rownames(pca_plot_data)), vjust = -1, color = "black") +
  theme_bw(base_size = 14) +
  scale_color_manual(values = c("Mutant" = "red", "Wild-Type" = "green")) +
  labs(
    title = " (PCA)",
    x = paste0("PC1 (", round(percent_var[1], 1), "% variance)"),
    y = paste0("PC2 (", round(percent_var[2], 1), "% variance)"),
    color = "Group"
  ) +
  theme(legend.position = "bottom")


# 5.6  (Heatmap of All Proteins/Genes)
#  pheatmap （） heatmap_all 
heatmap_all <- pheatmap(df2,
         scale = 'row',
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         show_rownames = FALSE,
         silent = FALSE) #  silent = TRUE pheatmap，

# ---  ---
# 
row_order_all <- heatmap_all$tree_row$order

# ， 'df' 
df_heatmap_ordered <- df[row_order_all, ]

# Excel
write_xlsx(df_heatmap_ordered, "heatmap_all_proteins_data_ordered.xlsx")

# ，
print(" heatmap_all_proteins_data_ordered.xlsx ")
# ---  ---


# 6.  (Specific Gene Set Analysis)
# 'fatty_acid_related' sheet 
dg <- read_xlsx('analysis.xlsx', sheet = 'fatty_acid_related')

#  intensity_cols ，
dg_intensity <- dg[, intensity_cols]

#  pheatmap  heatmap_fatty_acid 
heatmap_fatty_acid <- pheatmap(dg_intensity,
         scale = 'row',
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         labels_row = dg$`Leading_gene`, # dg'Leading_gene'
         silent = FALSE) #  silent = TRUE pheatmap

# ---  ---
# 
row_order_fatty_acid <- heatmap_fatty_acid$tree_row$order

# ， 'dg' 
dg_heatmap_ordered <- dg[row_order_fatty_acid, ]

# Excel
write_xlsx(dg_heatmap_ordered, "heatmap_fatty_acid_data_ordered.xlsx")

# ，
print(" heatmap_fatty_acid_data_ordered.xlsx ")
# ---  ---

