import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import re

# --- 1. 数据加载与准备 ---
try:
    df = pd.read_csv('violin_plot.csv')
    print("✔️ 文件 'violin_plot.csv' 加载成功!")
except FileNotFoundError:
    print("❌ 错误：未找到 'violin_plot.csv' 文件。请确保它和您的代码在同一个文件夹中。")
    exit()

metabolite_col_name = df.columns[0]
df = df.set_index(metabolite_col_name)

# --- 2. 数据清洗 ---
print("\n🔍 正在检查并清洗数据...")
for col in df.columns:
    df[col] = pd.to_numeric(df[col], errors='coerce')
df.dropna(how='all', inplace=True)
print("✔️ 数据清洗完毕。")

# --- 3. 计算RSD的核心逻辑 ---
print("\n🔬 正在为每个代谢物计算RSD...")
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
print("✔️ RSD计算完成。")

qc_rsd_median = df_rsd[df_rsd['Category'] == 'QC']['RSD'].median()
sample_rsd_median = df_rsd[df_rsd['Category'] == 'Sample']['RSD'].median()
print(f"\n📊 QC组RSD中位数: {qc_rsd_median:.2f}%")
print(f"   样本组RSD中位数: {sample_rsd_median:.2f}%")

# --- 4. 绘图 (修改后的核心部分) ---
plt.style.use('seaborn-v0_8-whitegrid')
fig, ax = plt.subplots(figsize=(9, 9)) # 稍微调整画布比例

# 步骤1: 绘制小提琴图，但移除内部的图 (inner=None)
# 我们让它半透明 (alpha=0.7) 以便看清后面的箱线图
sns.violinplot(
    data=df_rsd,
    x='Category',
    y='RSD',
    order=['QC', 'Sample'],
    palette={'QC': '#4682B4', 'Sample': '#FFA07A'},
    inner=None, # <-- 关键改动1：移除内部绘图
    alpha=0.7,  # 使小提琴半透明
    ax=ax
)

# 步骤2: 在小提琴图上层叠加一个清晰的箱线图
sns.boxplot(
    data=df_rsd,
    x='Category',
    y='RSD',
    order=['QC', 'Sample'],
    width=0.2,          # <-- 关键改动2：设置一个较窄的宽度
    boxprops={'facecolor':'white', 'edgecolor':'black', 'alpha':0.9}, # 设置箱体样式
    whiskerprops={'color':'black', 'linewidth':1.5},      # 设置须线样式
    capprops={'color':'black', 'linewidth':1.5},          # 设置顶盖线样式
    medianprops={'color':'red', 'linewidth':2},           # 突出中位数
    ax=ax
)

# --- 5. 添加注释和美化 ---
ax.axhline(y=20, color='grey', linestyle='--', linewidth=1.2, zorder=0)
ax.text(ax.get_xlim()[1], 20, ' 20% Threshold', va='center', ha='left', color='grey')
ax.axhline(y=30, color='darkgrey', linestyle='--', linewidth=1.2, zorder=0)
ax.text(ax.get_xlim()[1], 30, ' 30% Threshold', va='center', ha='left', color='darkgrey')

ax.set_title('Distribution of Metabolite RSDs', fontsize=18, fontweight='bold', pad=20)
ax.set_xlabel('Group', fontsize=14, labelpad=15)
ax.set_ylabel('RSD (%)', fontsize=14, labelpad=15)
ax.tick_params(axis='both', which='major', labelsize=12)

# 关键改动3: 移除了 set_ylim，让Y轴自动缩放以包含所有数据
# ax.set_ylim(0, ...) # 此行已被删除

plt.tight_layout()
plt.show()