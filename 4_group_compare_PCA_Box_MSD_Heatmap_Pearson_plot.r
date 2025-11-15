# 1. 加载程序包 (Load Libraries)
# 确保安装了所有需要的包
# install.packages(c("ggplot2", "readxl", "psych", "pheatmap", "limma", "reshape2", "writexl", "dplyr", "RColorBrewer"))
library(ggplot2)
library(readxl)
library(psych)
library(pheatmap)
library(limma)
library(reshape2) # 显式加载 melt 函数所在的包
library(writexl)  # 用于写入Excel文件
library(dplyr)    # 用于数据整理
library(RColorBrewer) # 用于生成更多颜色

# 2. 定义列名和参数 (Define Column Names and Parameters)
# --- 主要修改区域 ---
# 定义四个实验组的样本列名
mut_cols <- c("Mut_1", "Mut_2", "Mut_3")
wt_cols  <- c("WT_1", "WT_2", "WT_3")
ctl_cols <- c("Con_1", "Con_2", "Con_3")
# 新增第四组的列名 (请根据您的数据修改)
ng_cols  <- c("NewGroup_1", "NewGroup_2", "NewGroup_3")

# 将所有组信息整合到一个列表中，方便管理
group_list <- list(
  "Mutant" = mut_cols,
  "Wild-Type" = wt_cols,
  "Control" = ctl_cols,
  "NewGroup" = ng_cols  # 新增组
)

# 合并所有强度相关的列名
intensity_cols <- unlist(group_list)

# 定义组名和颜色
group_names <- names(group_list)
# 使用RColorBrewer来获取一组区分明显的颜色，也可以手动指定
# colors <- c("red", "green", "blue", "purple")
colors <- brewer.pal(length(group_names), "Set1")
names(colors) <- group_names


# 3. 数据读取与预处理 (Data Loading and Pre-processing)
# 假设您的 'lfq_analysis' sheet 包含了所有12个样本的强度数据
# 请确保Excel文件中有新增的第四组数据
df <- read_xlsx('analysis.xlsx', sheet = 'lfq_analysis')

# 提取强度数据用于后续所有QC图
# 这是一个关键步骤，后续所有图都将基于这个df2数据框
# 确保所有在 intensity_cols 中定义的列都存在于 df 中
if (!all(intensity_cols %in% colnames(df))) {
  stop("错误: Excel文件中缺少部分样本列，请检查 'intensity_cols' 和文件内容。")
}
df2 <- df[, intensity_cols]


# 4. 数据可视化 (Data Visualization)

# 4.1 样本强度分布箱线图 (Boxplot of Sample Intensities)
df2m <- reshape2::melt(df2, measure.vars = intensity_cols)
# 动态生成颜色映射
color_map_boxplot <- unlist(sapply(group_names, function(g) rep(colors[g], length(group_list[[g]]))))
names(color_map_boxplot) <- intensity_cols

ggplot(df2m, aes(x = variable, y = value, fill = variable)) +
  geom_boxplot() +
  scale_fill_manual(values = color_map_boxplot) +
  theme_bw() +
  theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "各样本强度分布箱线图", x = 'Sample', y = 'Log2 LFQ Intensity')

# 4.2 样本间相关性散点图 (Pairs Panels)
pairs.panels(df2)

# 4.3 多维尺度分析图 (MDS Plot)
# 动态生成MDS颜色
group_membership <- sapply(colnames(df2), function(col) {
  names(which(sapply(group_list, function(g) col %in% g)))
})
mds_colors <- colors[group_membership]

plotMDS(df2, col = mds_colors, labels = colnames(df2))
legend("topleft", legend = group_names, fill = colors[group_names])

# 4.4 主成分分析图 (PCA Plot)
pca_data <- t(df2)
pca_result <- prcomp(pca_data, scale. = TRUE)
pca_plot_data <- as.data.frame(pca_result$x)

# 动态生成分组因子
pca_plot_data$group <- factor(group_membership, levels = group_names)

percent_var <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100
ggplot(pca_plot_data, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 5, alpha = 0.8) +
  geom_text(aes(label = rownames(pca_plot_data)), vjust = -1.5, color = "black", size = 3) +
  theme_bw(base_size = 14) +
  scale_color_manual(values = colors) +
  labs(
    title = "主成分分析 (PCA)",
    x = paste0("PC1 (", round(percent_var[1], 1), "% variance)"),
    y = paste0("PC2 (", round(percent_var[2], 1), "% variance)"),
    color = "Group"
  ) +
  theme(legend.position = "bottom")

# 准备热图注释信息 (在两种热图前准备好)
group_vector <- factor(group_membership, levels = group_names)
annotation_col <- data.frame(Group = group_vector)
rownames(annotation_col) <- intensity_cols
annotation_colors <- list(Group = colors)


# 4.5 整体表达热图 (按表达谱聚类)
# 捕获 pheatmap 的输出对象
heatmap_all_clustered <- pheatmap(df2,
         scale = 'row',
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         cluster_cols = TRUE,          # <-- 启用列聚类
         cluster_rows = TRUE,
         show_rownames = FALSE,
         main = "全部蛋白/基因表达热图 (按表达谱聚类)",
         silent = TRUE) # 使用 silent=TRUE 避免在控制台重复绘图

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
         cluster_cols = FALSE,         # <-- 禁用列聚类
         cluster_rows = TRUE,
         show_rownames = FALSE,
         main = "全部蛋白/基因表达热图 (按分组排序)")


# 5. 特定基因集分析 (Specific Gene Set Analysis)
# 假设 'fatty_acid_related' 是您感兴趣的特定基因列表
dg <- read_xlsx('analysis.xlsx', sheet = 'fatty_acid_related')

# 检查特定基因集的数据框是否也包含所有样本列
if (all(intensity_cols %in% colnames(dg))) {
  dg_intensity <- dg[, intensity_cols]

  # 5.1 特定基因集热图 (按表达谱聚类)
  # 捕获 pheatmap 的输出对象
  heatmap_specific_clustered <- pheatmap(dg_intensity,
           scale = 'row',
           annotation_col = annotation_col,
           annotation_colors = annotation_colors,
           cluster_cols = TRUE,          # <-- 启用列聚类
           cluster_rows = TRUE,
           labels_row = dg$`Leading_gene`, # 确保您的dg数据框中有这一列
           main = "特定基因集热图 (按表达谱聚类)",
           silent = TRUE) # 使用 silent=TRUE 避免重复绘图

  # 从捕获的对象中提取行的顺序
  row_order_specific <- heatmap_specific_clustered$tree_row$order
  # 使用此顺序来排列原始的 'dg' 数据框
  dg_heatmap_ordered <- dg[row_order_specific, ]
  # 将排序后的数据框写入一个新的Excel文件
  write_xlsx(dg_heatmap_ordered, "heatmap_specific_genes_ordered.xlsx")
  # 打印消息，告知用户文件已保存
  print("与“特定基因集”热图行顺序对应的数据已保存到 heatmap_specific_genes_ordered.xlsx 文件中。")

  # 5.2 特定基因集热图 (按预设分组排序)
  pheatmap(dg_intensity,
           scale = 'row',
           annotation_col = annotation_col,
           annotation_colors = annotation_colors,
           cluster_cols = FALSE,         # <-- 禁用列聚类
           cluster_rows = TRUE,
           labels_row = dg$`Leading_gene`,
           main = "特定基因集表达热图 (按分组排序)")

} else {
  print("警告: 'fatty_acid_related' 工作表中缺少部分样本列，跳过特定基因集热图分析。")
}
