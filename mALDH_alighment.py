import os
import subprocess
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import AlignIO

# --- Step 1: Perform Multiple Sequence Alignment ---
# 定义输入和输出文件名
fasta_file = "uniprotkb_mALDH_sequence.fasta"
aligned_file = "aligned_sequences.aln"
clustal_format = "clustal"

print(f"Running Multiple Sequence Alignment on '{fasta_file}'...")

# 构建并执行Clustal Omega命令
# 注意: 你必须先在你的系统中安装 clustalo
# --force 覆盖已存在的文件, --outfmt 指定输出格式
try:
    subprocess.run(
        ["clustalo", "-i", fasta_file, "-o", aligned_file, "--outfmt=" + clustal_format, "--force"],
        check=True,
        capture_output=True, 
        text=True
    )
    print(f"Alignment complete. Output saved to '{aligned_file}'")
except FileNotFoundError:
    print("\n--- ERROR ---")
    print("Clustal Omega (clustalo) command not found.")
    print("Please install it on your system or use a web server like https://www.ebi.ac.uk/Tools/msa/clustalo/ to align your sequences.")
    print(f"Then, save the alignment in Clustal format as '{aligned_file}' and re-run this script.")
    exit()
except subprocess.CalledProcessError as e:
    print("\n--- ERROR during alignment ---")
    print(e.stderr)
    exit()


# --- Step 2: Read Alignment and Calculate Identity Matrix ---
print("Calculating identity matrix...")

# 读取比对结果
alignment = AlignIO.read(aligned_file, clustal_format)

# 提取简化的蛋白质名称用于标签
# 例如从 "sp|P24549|AL1A1_MOUSE" 提取 "AL1A1_MOUSE"
ids = [rec.id.split('|')[2] for rec in alignment]
num_sequences = len(ids)

# 初始化矩阵
identity_matrix = np.zeros((num_sequences, num_sequences))

# 逐对计算一致性
for i in range(num_sequences):
    for j in range(i, num_sequences):
        seq1 = alignment[i].seq
        seq2 = alignment[j].seq
        
        matches = 0
        aligned_positions = 0
        
        # 只比较两个序列中都不是gap的位置
        for k in range(alignment.get_alignment_length()):
            if seq1[k] != '-' and seq2[k] != '-':
                aligned_positions += 1
                if seq1[k] == seq2[k]:
                    matches += 1
        
        # 计算百分比
        if aligned_positions == 0:
            identity = 0
        else:
            identity = (matches / aligned_positions) * 100
        
        identity_matrix[i, j] = identity
        identity_matrix[j, i] = identity

# --- Step 3: Visualize the Heatmap ---
print("Generating heatmap...")

# 将矩阵转换为Pandas DataFrame以便添加标签
df_heatmap = pd.DataFrame(identity_matrix, index=ids, columns=ids)

# 设置绘图尺寸
plt.figure(figsize=(12, 10))

# 使用Seaborn绘制热图
sns.heatmap(
    df_heatmap, 
    cmap="Blues",       # 使用蓝色系，与你的示例图一致
    annot=False,        # 不在格子上显示数值
    vmin=10,            # 设置颜色条的最小值
    vmax=100            # 设置颜色条的最大值
)

# 添加标题和调整标签
plt.title("Protein Sequence Identity - Mouse ALDH Family", fontsize=16)
plt.xticks(rotation=90)
plt.yticks(rotation=0)
plt.tight_layout() # 自动调整布局

# --- 修改之处：将输出格式改为 SVG ---
# 定义SVG格式的输出文件名
output_image_svg = "mouse_aldh_identity_heatmap.svg"
# --- 备用方案：保存为 PDF ---
output_image_pdf = "mouse_aldh_identity_heatmap.pdf"
plt.savefig(output_image_pdf, format='pdf', bbox_inches='tight')

print(f"\nHeatmap successfully saved as '{output_image_pdf}'")
# 使用savefig函数保存为SVG格式
# Matplotlib会根据文件扩展名自动选择格式
plt.savefig(output_image_svg, format='svg')

print(f"Heatmap successfully saved as '{output_image_svg}'")