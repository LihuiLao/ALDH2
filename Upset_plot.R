# --- Upset Plot 脚本 ---

# 步骤 1: 安装和加载必要的包
# 检查 UpSetR 包是否已安装，如果没有则自动安装
if (!requireNamespace("UpSetR", quietly = TRUE)) {
  install.packages("UpSetR")
}
# 加载 UpSetR 包
library(UpSetR)

# 步骤 2: 读取你的数据文件
# 重要提示: 请确保 "1.xlsx - Sheet1.csv" 这个文件和你的R脚本在同一个文件夹下
# 或者，你也可以在这里使用文件的完整路径，例如 "C:/Users/12539/OneDrive/Desktop/1.xlsx - Sheet1.csv"
# R语言中的路径请使用正斜杠 "/"
file_path <- "1.csv"
metabolite_data <- read.csv(file_path, header = TRUE, sep = ",")

# 步骤 3: 准备数据以适应 UpSetR 格式
# UpSetR 的 fromExpression 函数需要一个只包含 0 和 1 的数据框
# 我们的数据中是 TRUE/FALSE，需要进行转换。同时，第一列是代谢物ID，不需要转换。
# 选取组别列（第2列到最后1列）
group_columns <- metabolite_data[, 2:ncol(metabolite_data)]

# 将 TRUE 转换为 1, FALSE 转换为 0
group_columns[group_columns == TRUE] <- 1
group_columns[group_columns == FALSE] <- 0

# 步骤 4: 绘制 Upset Plot
# 使用 upset() 函数创建图表
# nsets: 你想要展示的组别数量，这里是3个
# nintersects: 你想在图上展示多少个交集，可以设大一点比如40，NA表示全部展示
# order.by = "freq": 按照交集的大小（包含的代谢物数量）从大到小排序
# text.scale: 调整字体大小，可以根据你的需要修改
# point.size: 调整矩阵中圆点的大小
# line.size: 调整矩阵中连线的粗细
upset_plot <- upset(
  as.data.frame(group_columns), 
  nsets = 3, 
  nintersects = NA, 
  order.by = "freq",
  mb.ratio = c(0.6, 0.4),
  text.scale = c(1.5, 1.5, 1.2, 1.2, 1.5, 1.2), # 标题, x轴, y轴, 交集大小标题, 交集大小刻度, 组名
  point.size = 3.5,
  line.size = 2
)

# 步骤 5: 显示和保存图表
# 在RStudio的Plots窗口中直接显示图表
print(upset_plot)

# 如果需要将图表保存为PDF文件 (可选)
# pdf("Metabolite_Upset_Plot.pdf", width = 10, height = 7)
# print(upset_plot)
# dev.off() # 关闭PDF设备，完成保存

# --- 脚本结束 ---