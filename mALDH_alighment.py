import os
import subprocess
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import AlignIO

# --- Step 1: Perform Multiple Sequence Alignment ---
# 
fasta_file = "uniprotkb_mALDH_sequence.fasta"
aligned_file = "aligned_sequences.aln"
clustal_format = "clustal"

print(f"Running Multiple Sequence Alignment on '{fasta_file}'...")

# Clustal Omega
# :  clustalo
# --force , --outfmt 
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

# 
alignment = AlignIO.read(aligned_file, clustal_format)

# 
#  "sp|P24549|AL1A1_MOUSE"  "AL1A1_MOUSE"
ids = [rec.id.split('|')[2] for rec in alignment]
num_sequences = len(ids)

# 
identity_matrix = np.zeros((num_sequences, num_sequences))

# 
for i in range(num_sequences):
    for j in range(i, num_sequences):
        seq1 = alignment[i].seq
        seq2 = alignment[j].seq
        
        matches = 0
        aligned_positions = 0
        
        # gap
        for k in range(alignment.get_alignment_length()):
            if seq1[k] != '-' and seq2[k] != '-':
                aligned_positions += 1
                if seq1[k] == seq2[k]:
                    matches += 1
        
        # 
        if aligned_positions == 0:
            identity = 0
        else:
            identity = (matches / aligned_positions) * 100
        
        identity_matrix[i, j] = identity
        identity_matrix[j, i] = identity

# --- Step 3: Visualize the Heatmap ---
print("Generating heatmap...")

# Pandas DataFrame
df_heatmap = pd.DataFrame(identity_matrix, index=ids, columns=ids)

# 
plt.figure(figsize=(12, 10))

# Seaborn
sns.heatmap(
    df_heatmap, 
    cmap="Blues",       # ，
    annot=False,        # 
    vmin=10,            # 
    vmax=100            # 
)

# 
plt.title("Protein Sequence Identity - Mouse ALDH Family", fontsize=16)
plt.xticks(rotation=90)
plt.yticks(rotation=0)
plt.tight_layout() # 

# --- ： SVG ---
# SVG
output_image_svg = "mouse_aldh_identity_heatmap.svg"
# --- ： PDF ---
output_image_pdf = "mouse_aldh_identity_heatmap.pdf"
plt.savefig(output_image_pdf, format='pdf', bbox_inches='tight')

print(f"\nHeatmap successfully saved as '{output_image_pdf}'")
# savefigSVG
# Matplotlib
plt.savefig(output_image_svg, format='svg')

print(f"Heatmap successfully saved as '{output_image_svg}'")