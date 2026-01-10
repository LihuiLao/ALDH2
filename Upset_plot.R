# --- Upset Plot  ---

#  1: 
#  UpSetR ，
if (!requireNamespace("UpSetR", quietly = TRUE)) {
  install.packages("UpSetR")
}
#  UpSetR 
library(UpSetR)

#  2: 
# :  "1.xlsx - Sheet1.csv" R
# ，， "C:/Users/12539/OneDrive/Desktop/1.xlsx - Sheet1.csv"
# R "/"
file_path <- "1.csv"
metabolite_data <- read.csv(file_path, header = TRUE, sep = ",")

#  3:  UpSetR 
# UpSetR  fromExpression  0  1 
#  TRUE/FALSE，，ID，
# （21）
group_columns <- metabolite_data[, 2:ncol(metabolite_data)]

#  TRUE  1, FALSE  0
group_columns[group_columns == TRUE] <- 1
group_columns[group_columns == FALSE] <- 0

#  4:  Upset Plot
#  upset() 
# nsets: ，3
# nintersects: ，40，NA
# order.by = "freq": （）
# text.scale: ，
# point.size: 
# line.size: 
upset_plot <- upset(
  as.data.frame(group_columns), 
  nsets = 3, 
  nintersects = NA, 
  order.by = "freq",
  mb.ratio = c(0.6, 0.4),
  text.scale = c(1.5, 1.5, 1.2, 1.2, 1.5, 1.2), # , x, y, , , 
  point.size = 3.5,
  line.size = 2
)

#  5: 
# RStudioPlots
print(upset_plot)

# PDF ()
# pdf("Metabolite_Upset_Plot.pdf", width = 10, height = 7)
# print(upset_plot)
# dev.off() # PDF，

# ---  ---