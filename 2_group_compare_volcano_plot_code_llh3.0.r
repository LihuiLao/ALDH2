# 1. 加载程序库 (Load Libraries)
# 确保安装了所有需要的包
# install.packages(c("ggplot2", "readxl", "psych", "pheatmap", "limma", "reshape2", "writexl", "dplyr"))

library(ggplot2)
library(readxl)
library(psych)
library(pheatmap)
library(limma)
library(reshape2) # 显式加载 melt 函数所在的包
library(writexl)  # 新增：用于写入Excel文件
library(dplyr)    # 新增：用于数据整理，如此处的列重排

# 2. 定义列名和参数 (Define Column Names and Parameters)
# 定义实验组和对照组的样本列名
mut_cols <- c("Mut_1", "Mut_2", "Mut_3")
wt_cols <- c("WT_1", "WT_2", "WT_3")
# 合并所有强度相关的列名
intensity_cols <- c(mut_cols, wt_cols)
# 定义差异倍数所在的列名
fc_col <- "FC"
# 定义p值所在的列名
pval_col <- "p-value"

# 火山图曲线参数
s0 <- 0.1      # Fudge factor，调节曲线弯曲程度
f <- 5         # 自由度 (degrees of freedom) 自由度 = 样品数-组别数+1
confidence <- 0.95 # 置信水平
ta <- qt(confidence, f) # 根据置信度和自由度计算t分布的临界值

# 3. 数据加载与预处理 (Data Loading and Pre-processing)
# 请确保 'analysis.xlsx' 文件在您的R工作目录下，或者提供完整路径
df <- read_xlsx('analysis.xlsx', sheet = 'lfq_analysis')

# 定义火山图曲线函数
data_calsmoothcurve <- function(x, ta, s0, f) {
  y <- ta * (1 + (s0 / ((abs(x) / ta) - s0)))
  y <- -log10(2 * (1 - pt(y, df = f)))
  return(y)
}

# 计算-log10(p-value)
df$lgpvalue <- -log10(df[[pval_col]])

# 计算火山图曲线的y轴阈值
df$fudge <- apply(df, 1, function(x) {
  data_calsmoothcurve(as.numeric(x[fc_col]), ta, s0, f)
})

# 判断显著性
# 初始化'sig'列，'0'代表不显著
df$sig <- '0'

# 使用逻辑掩码（mask）进行判断
# 条件1: p-value足够小 (lgpvalue > fudge)
# 条件2: FC变化足够大 (FC > ta * s0 或 FC < -ta * s0)

# 标记上调的蛋白/基因 ('1')
up_regulated_mask <- (df$lgpvalue > df$fudge) & (df[[fc_col]] > ta * s0)
df$sig[up_regulated_mask] <- '1'

# 标记下调的蛋白/基因 ('-1')
down_regulated_mask <- (df$lgpvalue > df$fudge) & (df[[fc_col]] < -ta * s0)
df$sig[down_regulated_mask] <- '-1'

# (可选)将sig列转换为因子类型，便于后续绘图和分析
df$sig <- factor(df$sig, levels = c('-1', '0', '1'))


# 4. 导出火山图数据到Excel (Export Volcano Plot Data to Excel)
# 这里假设您的数据中有 'Protein IDs' 和 'Gene names' 列，请根据实际情况修改
# 如果没有这些列，可以注释掉 select() 这一行
df_to_export <- df %>%
  select(any_of(c("Protein IDs", "Gene names")), !!fc_col, !!pval_col, lgpvalue, sig, fudge, everything())

# 将包含所有计算列（lgpvalue, fudge, sig）的数据框写入一个新的Excel文件
# 这个文件会保存在你的R工作目录下
write_xlsx(df_to_export, "volcano_analysis_output.xlsx")

# 打印消息，告知用户文件已保存
print("火山图分析的完整数据已保存到 volcano_analysis_output.xlsx 文件中。")


# 5. 数据可视化 (Data Visualization)

# 5.1 火山图 (Volcano Plot)
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


# 5.2 样本强度箱线图 (Boxplot of Sample Intensities)
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


# 5.3 样本相关性散点图矩阵 (Pairs Panels for Sample Correlation)
pairs.panels(df2)


# 5.4 多维尺度分析图 (MDS Plot)
mds_colors <- ifelse(colnames(df2) %in% mut_cols, "red", "green")
plotMDS(df2, col = mds_colors, labels = colnames(df2))
legend("topleft", legend = c("Mutant", "Wild-Type"), fill = c("red", "green"))

# 5.5 主成分分析图 (PCA Plot)
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
    title = "主成分分析 (PCA)",
    x = paste0("PC1 (", round(percent_var[1], 1), "% variance)"),
    y = paste0("PC2 (", round(percent_var[2], 1), "% variance)"),
    color = "Group"
  ) +
  theme(legend.position = "bottom")


# 5.6 整体表达热图 (Heatmap of All Proteins/Genes)
# 运行 pheatmap 并将结果（包含聚类信息）保存到变量 heatmap_all 中
heatmap_all <- pheatmap(df2,
         scale = 'row',
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         show_rownames = FALSE,
         silent = FALSE) # 使用 silent = TRUE 来阻止pheatmap直接打印图形，如果需要的话

# --- 新增代码开始 ---
# 从保存的热图对象中提取行的顺序
row_order_all <- heatmap_all$tree_row$order

# 根据热图的行顺序，对原始的完整数据框 'df' 进行排序
df_heatmap_ordered <- df[row_order_all, ]

# 将排序后的数据框写入一个新的Excel文件
write_xlsx(df_heatmap_ordered, "heatmap_all_proteins_data_ordered.xlsx")

# 打印消息，告知用户文件已保存
print("与整体表达热图行顺序对应的数据已保存到 heatmap_all_proteins_data_ordered.xlsx 文件中。")
# --- 新增代码结束 ---


# 6. 特定基因集分析 (Specific Gene Set Analysis)
# 'fatty_acid_related' sheet 也有相同的强度列
dg <- read_xlsx('analysis.xlsx', sheet = 'fatty_acid_related')

# 使用 intensity_cols 选择列，确保列顺序变化也能正确选择
dg_intensity <- dg[, intensity_cols]

# 运行 pheatmap 并将结果保存到变量 heatmap_fatty_acid 中
heatmap_fatty_acid <- pheatmap(dg_intensity,
         scale = 'row',
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         labels_row = dg$`Leading_gene`, # 假设dg文件中有'Leading_gene'列
         silent = FALSE) # 使用 silent = TRUE 来阻止pheatmap直接打印图形

# --- 新增代码开始 ---
# 从保存的热图对象中提取行的顺序
row_order_fatty_acid <- heatmap_fatty_acid$tree_row$order

# 根据热图的行顺序，对原始的 'dg' 数据框进行排序
dg_heatmap_ordered <- dg[row_order_fatty_acid, ]

# 将排序后的数据框写入一个新的Excel文件
write_xlsx(dg_heatmap_ordered, "heatmap_fatty_acid_data_ordered.xlsx")

# 打印消息，告知用户文件已保存
print("与脂肪酸相关蛋白热图行顺序对应的数据已保存到 heatmap_fatty_acid_data_ordered.xlsx 文件中。")
# --- 新增代码结束 ---

