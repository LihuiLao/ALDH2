# 加载必要的包
library(ggplot2)      # 必须加载！ggtitle函数在这个包里
library(ggseqlogo)
library(Biostrings)
# 读取FASTA文件
seqs <- readAAStringSet("C:/Users/12539/OneDrive/Desktop/alighment.fasta")
seqs_char <- as.character(seqs)
# 生成logo图
p <- ggseqlogo(seqs_char, method='bits', seq_type='aa') +
  ggtitle('hALDH Amino Acid Conservation') +
  theme_logo() +
  theme(axis.text.x = element_text(size=10))
# 保存为SVG
ggsave("C:/Users/12539/OneDrive/Desktop/hALDH_logo.svg", 
       p, width=15, height=5, units="in")