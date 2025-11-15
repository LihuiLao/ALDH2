import pandas as pd
import numpy as np
import time

def extract_peaks_with_conditions_optimized_corrected(file_path, mass_difference=2.00671, mz_tolerance=0.02, retention_tolerance=0.05, intensity_ratio_range=(1, 2.0)):
    """
    从Excel文件中提取符合特定条件的峰对（优化并修正版）。
    修正点：
    - np.searchsorted 的边界条件确保与 < mz_tolerance 一致。
    - 平均强度NumPy数组在数据排序后创建。
    """
    print("开始读取Excel文件...")
    data = pd.read_excel(file_path)
    print(f"成功读取数据，共 {len(data)} 行。")

    # --- 定义列名 ---
    peak_name_col = 'Peak Name'
    mz_col = 'm/z'
    rt_col = 'Retention time'
    
    # --- 检查必需列 ---
    original_intensity_cols_required = []
    for i in range(1, 28):
        for j in range(1, 4):
            original_intensity_cols_required.append(f'Intensity_LLH_{i}_{j}')
    required_columns_check = [peak_name_col, mz_col, rt_col] + original_intensity_cols_required
    missing_cols = [col for col in required_columns_check if col not in data.columns]
    if missing_cols:
        raise ValueError(f"Excel表格必须包含以下列: '{peak_name_col}', '{mz_col}', '{rt_col}' "
                         f"以及所有形如 'Intensity_LLH_X_Y' 的强度列. 缺失的列: {missing_cols}")

    print("正在计算每个样本的平均强度并添加到DataFrame...")
    mean_intensity_col_names_list = [] # 用于存储平均强度列的名称
    for i in range(1, 28):
        replicate_cols = [f'Intensity_LLH_{i}_{j}' for j in range(1, 4)]
        mean_col_name = f'Mean_Intensity_LLH_{i}'
        mean_intensity_col_names_list.append(mean_col_name)
        if all(col in data.columns for col in replicate_cols):
            data[mean_col_name] = data[replicate_cols].mean(axis=1)
        else:
            raise ValueError(f"计算平均强度时，样本 {i} 的部分重复列缺失: {replicate_cols}")
    print(f"成功计算 {len(mean_intensity_col_names_list)} 个平均强度列。")

    print("正在排序数据...")
    data = data.sort_values(by=mz_col).reset_index(drop=True)

    print("正在将排好序的列转换为NumPy数组...")
    # 将核心列转换为NumPy数组以提高访问速度 (在排序后进行)
    mz_array = data[mz_col].to_numpy()
    rt_array = data[rt_col].to_numpy()
    peak_name_array = data[peak_name_col].to_numpy()

    # 存储平均强度列的NumPy数组 (在排序后进行)
    mean_intensity_col_data_np = {} 
    for col_name in mean_intensity_col_names_list:
        mean_intensity_col_data_np[col_name] = data[col_name].to_numpy()

    # 存储原始强度列的NumPy数组 (在排序后进行)
    original_intensity_arrays_np = {}
    for s_idx in range(1, 28):
        for r_idx in range(1, 4):
            col_name = f'Intensity_LLH_{s_idx}_{r_idx}'
            if col_name in data.columns: 
                 original_intensity_arrays_np[col_name] = data[col_name].to_numpy()

    extracted_peaks = []
    num_rows = len(data)
    
    if num_rows < 2:
        print("数据行数不足2行，无法进行峰对比较。")
        return pd.DataFrame()
        
    actual_comparisons_made = 0 
    start_time = time.time()

    print("开始使用优化方法寻找符合条件的峰对...")
    for i in range(num_rows - 1): 
        mz_i = mz_array[i]
        rt_i = rt_array[i]

        # 定义行j的m/z搜索窗口，基于行i的m/z和mass_difference
        # 我们需要 L < mz_j < U
        # L = mz_i + mass_difference - mz_tolerance
        # U = mz_i + mass_difference + mz_tolerance
        
        # 修正点1: searchsorted的side参数调整以匹配严格小于的逻辑
        lower_mz_for_j_exclusive = mz_i + mass_difference - mz_tolerance
        upper_mz_for_j_exclusive = mz_i + mass_difference + mz_tolerance

        # 找到第一个 mz_array[j_idx] > lower_mz_for_j_exclusive 的索引
        start_j_candidate_idx = np.searchsorted(mz_array, lower_mz_for_j_exclusive, side='right')
        # 找到第一个 mz_array[j_idx] >= upper_mz_for_j_exclusive 的索引
        end_j_candidate_idx = np.searchsorted(mz_array, upper_mz_for_j_exclusive, side='left')
        
        for j_idx in range(start_j_candidate_idx, end_j_candidate_idx):
            if j_idx <= i: 
                continue
            
            actual_comparisons_made += 1

            # m/z 条件通过 searchsorted 的窗口选择已满足 (L < mz_j < U)
            # 现在检查保留时间差异
            rt_j = rt_array[j_idx]
            retention_diff = abs(rt_j - rt_i)

            if retention_diff < retention_tolerance:
                at_least_one_ratio_in_range = False
                for k_sample_num in range(1, 28): 
                    mean_intensity_col_name = f'Mean_Intensity_LLH_{k_sample_num}'
                    
                    # 修正点2: 从排序后的NumPy数组中获取平均强度
                    intensity1 = mean_intensity_col_data_np[mean_intensity_col_name][i]
                    intensity2 = mean_intensity_col_data_np[mean_intensity_col_name][j_idx]
                    
                    if pd.notna(intensity1) and pd.notna(intensity2) and intensity1 > 0 and intensity2 > 0:
                        ratio = max(intensity1, intensity2) / min(intensity1, intensity2)
                        if intensity_ratio_range[0] <= ratio <= intensity_ratio_range[1]:
                            at_least_one_ratio_in_range = True
                            break 
                
                if at_least_one_ratio_in_range:
                    peak_pair = {
                        'Peak 1 Peak Name': peak_name_array[i],
                        'Peak 1 m/z': mz_i,
                        'Peak 1 Retention time': rt_i,
                        'Peak 2 Peak Name': peak_name_array[j_idx],
                        'Peak 2 m/z': mz_array[j_idx],
                        'Peak 2 Retention time': rt_j,
                        'm/z Difference': abs(mz_array[j_idx] - mz_i), 
                        'Retention Time Difference': retention_diff
                    }
                    
                    for s_idx in range(1, 28):
                        for r_idx in range(1, 4):
                            orig_col_name = f'Intensity_LLH_{s_idx}_{r_idx}'
                            peak_pair[f'Peak 1 {orig_col_name}'] = original_intensity_arrays_np[orig_col_name][i]
                            peak_pair[f'Peak 2 {orig_col_name}'] = original_intensity_arrays_np[orig_col_name][j_idx]
                        
                        mean_col = f'Mean_Intensity_LLH_{s_idx}'
                        peak_pair[f'Peak 1 {mean_col}'] = mean_intensity_col_data_np[mean_col][i]
                        peak_pair[f'Peak 2 {mean_col}'] = mean_intensity_col_data_np[mean_col][j_idx]
                        
                    extracted_peaks.append(peak_pair)
        
        if (i + 1) % 1000 == 0 or i == num_rows - 2 : 
            elapsed_time = time.time() - start_time
            progress_outer = (i + 1) / num_rows * 100
            print(f"外层循环进度: {progress_outer:.2f}% ({i+1}/{num_rows}), "
                  f"已进行 {actual_comparisons_made} 次实际峰对比较, 用时: {elapsed_time:.2f} 秒")

    print(f"峰对比较完成，共找到 {len(extracted_peaks)} 对符合所有条件的峰。")
    print(f"总共进行了 {actual_comparisons_made} 次实际的峰对比较。")
    if not extracted_peaks:
        print("未找到符合条件的峰对。")
        return pd.DataFrame() 
        
    return pd.DataFrame(extracted_peaks)

# --- 使用示例 ---
file_path = 'C:/Users/12539/OneDrive/Desktop/OH_Height_extraction.xlsx' 

print("脚本开始执行 (优化并修正版)...")
try:
    result_df = extract_peaks_with_conditions_optimized_corrected(
        file_path,
        mass_difference=2.00671,
        mz_tolerance=0.02,
        retention_tolerance=0.05,
        intensity_ratio_range=(1, 2.0) 
    )

    if not result_df.empty:
        output_file = 'OH_extracted_peaks.xlsx' 
        result_df.to_excel(output_file, index=False)
        print(f"脚本执行完毕，结果已保存到 {output_file}")
    else:
        print("脚本执行完毕，但没有找到符合条件的峰对，未生成结果文件。")

except FileNotFoundError:
    print(f"错误: 文件未找到 {file_path}。请检查路径是否正确。")
except ValueError as ve:
    print(f"数据错误或配置错误: {ve}")
except Exception as e:
    print(f"发生未知错误: {e}")
    import traceback
    print("详细错误追踪:")
    traceback.print_exc()

