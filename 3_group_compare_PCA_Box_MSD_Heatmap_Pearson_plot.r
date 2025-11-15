# 1. 加载程序包 (Load Libraries)
# 确保安装了所有需要的包
# install.packages(c("ggplot2", "readxl", "psych", "pheatmap", "limma", "reshape2", "writexl", "dplyr"))
library(ggplot2)
library(readxl)
library(psych)
library(pheatmap)
library(limma)
library(reshape2) # 显式加载 melt 函数所在的包
library(writexl)  # 新增：用于写入Excel文件
library(dplyr)    # 新增：用于数据整理

# 2. 定义列名和参数 (Define Column Names and Parameters)
# 定义三个实验组的样本列名
mut_cols <- c("Mut_1", "Mut_2", "Mut_3")
wt_cols <- c("WT_1", "WT_2", "WT_3")
ctl_cols <- c("Con_1", "Con_2", "Con_3")

# 合并所有强度相关的列名
intensity_cols <- c(mut_cols, wt_cols, ctl_cols)

# 正确的多组比较应使用ANOVA等方法，这里我们专注于探索性数据分析图。

# 3. 数据读取与预处理 (Data Loading and Pre-processing)
# 假设你的 'lfq_analysis' sheet 包含了所有9个样本的强度数据
df <- read_xlsx('analysis.xlsx', sheet = 'lfq_analysis')

# 提取强度数据用于后续所有QC图
# 这是一个关键步骤，后续所有图都将基于这个df2数据框
df2 <- df[, intensity_cols]

# 4. 数据可视化 (Data Visualization)

# 4.1 样本强度分布箱线图 (Boxplot of Sample Intensities)
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
  labs(title = "各样本强度分布箱线图", x = 'Sample', y = 'Log2 LFQ Intensity')

# 4.2 样本间相关性散点图 (Pairs Panels)
pairs.panels(df2)

# 4.3 多维尺度分析图 (MDS Plot)
mds_colors <- ifelse(colnames(df2) %in% mut_cols, "red",
                      ifelse(colnames(df2) %in% wt_cols, "green", "blue"))
plotMDS(df2, col = mds_colors, labels = colnames(df2))
legend("topleft", legend = c("Mutant", "Wild-Type", "Control"), fill = c("red", "green", "blue"))

# 4.4 主成分分析图 (PCA Plot)
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
    title = "主成分分析 (PCA)",
    x = paste0("PC1 (", round(percent_var[1], 1), "% variance)"),
    y = paste0("PC2 (", round(percent_var[2], 1), "% variance)"),
    color = "Group"
  ) +
  theme(legend.position = "bottom")

# 准备热图注释信息 (在两种热图前准备好)
# 确保因子水平顺序与 intensity_cols 一致
group_vector <- rep(c("Mutant", "Wild-Type", "Control"), c(length(mut_cols), length(wt_cols), length(ctl_cols)))
names(group_vector) <- c(mut_cols, wt_cols, ctl_cols)
group_vector <- group_vector[intensity_cols] # 按照 intensity_cols 的顺序排序

annotation_col <- data.frame(Group = factor(group_vector, levels = c("Mutant", "Wild-Type", "Control")))
rownames(annotation_col) <- intensity_cols
annotation_colors <- list(Group = c(Mutant = "red", `Wild-Type` = "green", Control = "blue"))


# 4.5 整体表达热图 (按表达谱聚类)
# 捕获 pheatmap 的输出对象
heatmap_all_clustered <- pheatmap(df2,
         scale = 'row',
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         cluster_cols = TRUE,               # <-- 启用列聚类
         cluster_rows = TRUE,
         show_rownames = FALSE,
         main = "全部蛋白/基因表达热图 (按表达谱聚类)",
         silent = FALSE) # 使用 silent=TRUE 避免重复绘图

# 从捕获的对象中提取行的顺序
row_order_all <- heatmap_all_clustered$tree_row$order
# 使用此顺序来排列原始的完整数据框 'df'
df_heatmap_ordered <- df[row_order_all, ]
# 将排序后的数据框写入一个新的Excel文件
write_xlsx(df_heatmap_ordered, "heatmap_all_proteins_ordered.xlsx")
# 打印消息，告知用户文件已保存
print("与“全部蛋白”热图行顺序对应的数据已保存到 heatmap_all_proteins_ordered.xlsx 文件中。")


# 4.6 整体表达热图 (按预设分组排序)
pheatmap(df2,
         scale = 'row',
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         cluster_cols = FALSE,              # <-- 禁用列聚类
         cluster_rows = TRUE,
         show_rownames = FALSE,
         main = "全部蛋白/基因表达热图 (按分组排序)")


# 5. 特定基因集分析 (Specific Gene Set Analysis)
dg <- read_xlsx('analysis.xlsx', sheet = 'fatty_acid_related')
if (all(intensity_cols %in% colnames(dg))) {
  dg_intensity <- dg[, intensity_cols]

  # 5.1 脂肪酸相关蛋白热图 (按表达谱聚类)
  # 捕获 pheatmap 的输出对象
  heatmap_fa_clustered <- pheatmap(dg_intensity,
           scale = 'row',
           annotation_col = annotation_col,
           annotation_colors = annotation_colors,
           cluster_cols = TRUE,               # <-- 启用列聚类
           cluster_rows = TRUE,
           labels_row = dg$`Leading_gene`,
           main = "脂肪酸相关蛋白热图 (按表达谱聚类)",
           silent = FALSE) # 使用 silent=TRUE 避免重复绘图

  # 从捕获的对象中提取行的顺序
  row_order_fa <- heatmap_fa_clustered$tree_row$order
  # 使用此顺序来排列原始的 'dg' 数据框
  dg_heatmap_ordered <- dg[row_order_fa, ]
  # 将排序后的数据框写入一个新的Excel文件
  write_xlsx(dg_heatmap_ordered, "heatmap_fatty_acid_ordered.xlsx")
  # 打印消息，告知用户文件已保存
  print("与“脂肪酸相关蛋白”热图行顺序对应的数据已保存到 heatmap_fatty_acid_ordered.xlsx 文件中。")

  # 5.2 脂肪酸相关蛋白热图 (按预设分组排序)
  pheatmap(dg_intensity,
           scale = 'row',
           annotation_col = annotation_col,
           annotation_colors = annotation_colors,
           cluster_cols = FALSE,              # <-- 禁用列聚类
           cluster_rows = TRUE,
           labels_row = dg$`Leading_gene`,
           main = "脂肪酸相关蛋白表达热图 (按分组排序)")

} else {
  print("警告: 'fatty_acid_related' 工作表中缺少部分样本列，跳过此热图。")
}
