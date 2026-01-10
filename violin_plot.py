import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import re

# --- 1.  ---
try:
    df = pd.read_csv('violin_plot.csv')
    print("✔️  'violin_plot.csv' !")
except FileNotFoundError:
    print("❌ ： 'violin_plot.csv' ")
    exit()

metabolite_col_name = df.columns[0]
df = df.set_index(metabolite_col_name)

# --- 2.  ---
print("\n🔍 ...")
for col in df.columns:
    df[col] = pd.to_numeric(df[col], errors='coerce')
df.dropna(how='all', inplace=True)
print("✔️ ")

# --- 3. RSD ---
print("\n🔬 RSD...")
groups = {}
for col in df.columns:
    match = re.search(r'(QC|LLH_\d+)', col)
    if match:
        group_name = match.group(1)
        if group_name not in groups:
            groups[group_name] = []
        groups[group_name].append(col)

rsd_results = []
for group_name, columns in groups.items():
    if len(columns) < 2:
        continue
    
    group_df_log = df[columns]
    group_df_linear = np.power(2, group_df_log)
    mean_linear = group_df_linear.mean(axis=1)
    std_linear = group_df_linear.std(axis=1)
    rsd = np.divide(std_linear, mean_linear, out=np.zeros_like(std_linear), where=(mean_linear != 0)) * 100
    
    temp_df = pd.DataFrame({
        'RSD': rsd.dropna(),
        'Category': 'QC' if 'QC' in group_name else 'Sample'
    })
    rsd_results.append(temp_df)

df_rsd = pd.concat(rsd_results, ignore_index=True)
print("✔️ RSD")

qc_rsd_median = df_rsd[df_rsd['Category'] == 'QC']['RSD'].median()
sample_rsd_median = df_rsd[df_rsd['Category'] == 'Sample']['RSD'].median()
print(f"\n📊 QCRSD: {qc_rsd_median:.2f}%")
print(f"   RSD: {sample_rsd_median:.2f}%")

# --- 4.  () ---
plt.style.use('seaborn-v0_8-whitegrid')
fig, ax = plt.subplots(figsize=(9, 9)) # 

# 1: ， (inner=None)
#  (alpha=0.7) 
sns.violinplot(
    data=df_rsd,
    x='Category',
    y='RSD',
    order=['QC', 'Sample'],
    palette={'QC': '#4682B4', 'Sample': '#FFA07A'},
    inner=None, # <-- 1：
    alpha=0.7,  # 
    ax=ax
)

# 2: 
sns.boxplot(
    data=df_rsd,
    x='Category',
    y='RSD',
    order=['QC', 'Sample'],
    width=0.2,          # <-- 2：
    boxprops={'facecolor':'white', 'edgecolor':'black', 'alpha':0.9}, # 
    whiskerprops={'color':'black', 'linewidth':1.5},      # 
    capprops={'color':'black', 'linewidth':1.5},          # 
    medianprops={'color':'red', 'linewidth':2},           # 
    ax=ax
)

# --- 5.  ---
ax.axhline(y=20, color='grey', linestyle='--', linewidth=1.2, zorder=0)
ax.text(ax.get_xlim()[1], 20, ' 20% Threshold', va='center', ha='left', color='grey')
ax.axhline(y=30, color='darkgrey', linestyle='--', linewidth=1.2, zorder=0)
ax.text(ax.get_xlim()[1], 30, ' 30% Threshold', va='center', ha='left', color='darkgrey')

ax.set_title('Distribution of Metabolite RSDs', fontsize=18, fontweight='bold', pad=20)
ax.set_xlabel('Group', fontsize=14, labelpad=15)
ax.set_ylabel('RSD (%)', fontsize=14, labelpad=15)
ax.tick_params(axis='both', which='major', labelsize=12)

# 3:  set_ylim，Y
# ax.set_ylim(0, ...) # 

plt.tight_layout()
plt.show()