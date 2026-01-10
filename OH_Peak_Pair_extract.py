import pandas as pd
import numpy as np
import time

def extract_peaks_with_conditions_optimized_corrected(file_path, mass_difference=2.00671, mz_tolerance=0.02, retention_tolerance=0.05, intensity_ratio_range=(1, 2.0)):
    """
    Excel（）
    ：
    - np.searchsorted  < mz_tolerance 
    - NumPy
    """
    print("Excel...")
    data = pd.read_excel(file_path)
    print(f"， {len(data)} ")

    # ---  ---
    peak_name_col = 'Peak Name'
    mz_col = 'm/z'
    rt_col = 'Retention time'
    
    # ---  ---
    original_intensity_cols_required = []
    for i in range(1, 28):
        for j in range(1, 4):
            original_intensity_cols_required.append(f'Intensity_LLH_{i}_{j}')
    required_columns_check = [peak_name_col, mz_col, rt_col] + original_intensity_cols_required
    missing_cols = [col for col in required_columns_check if col not in data.columns]
    if missing_cols:
        raise ValueError(f"Excel: '{peak_name_col}', '{mz_col}', '{rt_col}' "
                         f" 'Intensity_LLH_X_Y' . : {missing_cols}")

    print("DataFrame...")
    mean_intensity_col_names_list = [] # 
    for i in range(1, 28):
        replicate_cols = [f'Intensity_LLH_{i}_{j}' for j in range(1, 4)]
        mean_col_name = f'Mean_Intensity_LLH_{i}'
        mean_intensity_col_names_list.append(mean_col_name)
        if all(col in data.columns for col in replicate_cols):
            data[mean_col_name] = data[replicate_cols].mean(axis=1)
        else:
            raise ValueError(f"， {i} : {replicate_cols}")
    print(f" {len(mean_intensity_col_names_list)} ")

    print("...")
    data = data.sort_values(by=mz_col).reset_index(drop=True)

    print("NumPy...")
    # NumPy ()
    mz_array = data[mz_col].to_numpy()
    rt_array = data[rt_col].to_numpy()
    peak_name_array = data[peak_name_col].to_numpy()

    # NumPy ()
    mean_intensity_col_data_np = {} 
    for col_name in mean_intensity_col_names_list:
        mean_intensity_col_data_np[col_name] = data[col_name].to_numpy()

    # NumPy ()
    original_intensity_arrays_np = {}
    for s_idx in range(1, 28):
        for r_idx in range(1, 4):
            col_name = f'Intensity_LLH_{s_idx}_{r_idx}'
            if col_name in data.columns: 
                 original_intensity_arrays_np[col_name] = data[col_name].to_numpy()

    extracted_peaks = []
    num_rows = len(data)
    
    if num_rows < 2:
        print("2，")
        return pd.DataFrame()
        
    actual_comparisons_made = 0 
    start_time = time.time()

    print("...")
    for i in range(num_rows - 1): 
        mz_i = mz_array[i]
        rt_i = rt_array[i]

        # jm/z，im/zmass_difference
        #  L < mz_j < U
        # L = mz_i + mass_difference - mz_tolerance
        # U = mz_i + mass_difference + mz_tolerance
        
        # 1: searchsortedside
        lower_mz_for_j_exclusive = mz_i + mass_difference - mz_tolerance
        upper_mz_for_j_exclusive = mz_i + mass_difference + mz_tolerance

        #  mz_array[j_idx] > lower_mz_for_j_exclusive 
        start_j_candidate_idx = np.searchsorted(mz_array, lower_mz_for_j_exclusive, side='right')
        #  mz_array[j_idx] >= upper_mz_for_j_exclusive 
        end_j_candidate_idx = np.searchsorted(mz_array, upper_mz_for_j_exclusive, side='left')
        
        for j_idx in range(start_j_candidate_idx, end_j_candidate_idx):
            if j_idx <= i: 
                continue
            
            actual_comparisons_made += 1

            # m/z  searchsorted  (L < mz_j < U)
            # 
            rt_j = rt_array[j_idx]
            retention_diff = abs(rt_j - rt_i)

            if retention_diff < retention_tolerance:
                at_least_one_ratio_in_range = False
                for k_sample_num in range(1, 28): 
                    mean_intensity_col_name = f'Mean_Intensity_LLH_{k_sample_num}'
                    
                    # 2: NumPy
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
            print(f": {progress_outer:.2f}% ({i+1}/{num_rows}), "
                  f" {actual_comparisons_made} , : {elapsed_time:.2f} ")

    print(f"， {len(extracted_peaks)} ")
    print(f" {actual_comparisons_made} ")
    if not extracted_peaks:
        print("")
        return pd.DataFrame() 
        
    return pd.DataFrame(extracted_peaks)

# ---  ---
file_path = 'C:/Users/12539/OneDrive/Desktop/OH_Height_extraction.xlsx' 

print(" ()...")
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
        print(f"， {output_file}")
    else:
        print("，，")

except FileNotFoundError:
    print(f":  {file_path}")
except ValueError as ve:
    print(f": {ve}")
except Exception as e:
    print(f": {e}")
    import traceback
    print(":")
    traceback.print_exc()

