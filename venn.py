import pandas as pd
from matplotlib_venn import venn3
import matplotlib.pyplot as plt

# CSV
df = pd.read_csv('Venn.csv')

# 1，ID
set1 = set(df[df['1_10_sig'] == 1]['Sample'])
set2 = set(df[df['4_10_sig'] == 1]['Sample'])
set3 = set(df[df['11_25_sig'] == 1]['Sample'])

# 
plt.figure(figsize=(10, 8))
venn3(
    [set1, set2, set3],
    set_labels=('Group 1_10_sig', 'Group 4_10_sig', 'Group 11_25_sig'),
    alpha=0.7
)

# 
plt.title("Venn Diagram of Sample Overlaps")
plt.savefig("Venn_Diagram.png", dpi=300) # 
plt.show() # 