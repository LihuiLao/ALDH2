# 
library(ggplot2)      # ！ggtitle
library(ggseqlogo)
library(Biostrings)
# FASTA
seqs <- readAAStringSet("C:/Users/12539/OneDrive/Desktop/alighment.fasta")
seqs_char <- as.character(seqs)
# logo
p <- ggseqlogo(seqs_char, method='bits', seq_type='aa') +
  ggtitle('hALDH Amino Acid Conservation') +
  theme_logo() +
  theme(axis.text.x = element_text(size=10))
# SVG
ggsave("C:/Users/12539/OneDrive/Desktop/hALDH_logo.svg", 
       p, width=15, height=5, units="in")