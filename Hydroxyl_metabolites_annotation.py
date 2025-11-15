import pandas as pd
import os
from datetime import datetime

def annotate_all_matches_advanced(raw_file, annotation_file, output_file, 
                                  mw_tolerance=0.03, 
                                  add_original_index=True,
                                  highlight_best_match=True):
    """
    增强版：添加更多功能
    
    新增参数:
        add_original_index: 是否添加原始行号
        highlight_best_match: 是否标记最佳匹配（MW最接近的）
    """
    
    print("=" * 60)
    print("Excel 文件匹配与注释工具 (增强版)")
    print("=" * 60)
    
    # 检查文件
    if not os.path.exists(raw_file):
        print(f"错误: 找不到文件 {raw_file}")
        return
    if not os.path.exists(annotation_file):
        print(f"错误: 找不到文件 {annotation_file}")
        return
    
    # 读取文件
    print(f"\n正在读取文件...")
    try:
        raw_df = pd.read_excel(raw_file)
        print(f"✓ 已读取 {raw_file}: {len(raw_df)} 行, {len(raw_df.columns)} 列")
        
        annotation_df = pd.read_excel(annotation_file)
        print(f"✓ 已读取 {annotation_file}: {len(annotation_df)} 行, {len(annotation_df.columns)} 列")
    except Exception as e:
        print(f"错误: 读取文件失败 - {e}")
        return
    
    # 检查MW列
    if 'MW' not in raw_df.columns or 'MW' not in annotation_df.columns:
        print("错误: 缺少 'MW' 列")
        return
    
    print(f"\n开始匹配处理...")
    print("-" * 60)
    
    results = []
    stats = {
        'total_rows': len(annotation_df),
        'no_match': 0,
        'single_match': 0,
        'multiple_match': 0,
    }
    
    # 遍历处理
    for idx, row in annotation_df.iterrows():
        mw_value = row['MW']
        
        # 查找匹配
        matches = raw_df[
            (raw_df['MW'] >= mw_value - mw_tolerance) & 
            (raw_df['MW'] <= mw_value + mw_tolerance)
        ].copy()
        
        # 计算MW差值并排序
        if not matches.empty:
            matches['abs_mw_diff'] = abs(matches['MW'] - mw_value)
            matches = matches.sort_values('abs_mw_diff')
        
        if matches.empty:
            # 无匹配
            stats['no_match'] += 1
            result_row = row.to_dict()
            if add_original_index:
                result_row['original_row_number'] = idx + 2  # Excel从第2行开始（有表头）
            result_row['match_status'] = 'No Match'
            result_row['match_index'] = None
            result_row['match_count'] = 0
            result_row['is_best_match'] = None
            results.append(result_row)
            
        elif len(matches) == 1:
            # 单个匹配
            stats['single_match'] += 1
            match_row = matches.iloc[0]
            result_row = row.to_dict()
            
            if add_original_index:
                result_row['original_row_number'] = idx + 2
            
            for col in raw_df.columns:
                if col != 'MW':
                    result_row[f'matched_{col}'] = match_row[col]
            
            result_row['matched_MW'] = match_row['MW']
            result_row['MW_difference'] = match_row['MW'] - mw_value
            result_row['abs_MW_difference'] = abs(match_row['MW'] - mw_value)
            result_row['match_status'] = 'Single Match'
            result_row['match_index'] = 1
            result_row['match_count'] = 1
            result_row['is_best_match'] = True
            
            results.append(result_row)
            
        else:
            # 多个匹配
            stats['multiple_match'] += 1
            
            for match_idx, (_, match_row) in enumerate(matches.iterrows(), 1):
                result_row = row.to_dict()
                
                if add_original_index:
                    result_row['original_row_number'] = idx + 2
                
                for col in raw_df.columns:
                    if col != 'MW':
                        result_row[f'matched_{col}'] = match_row[col]
                
                result_row['matched_MW'] = match_row['MW']
                result_row['MW_difference'] = match_row['MW'] - mw_value
                result_row['abs_MW_difference'] = abs(match_row['MW'] - mw_value)
                result_row['match_status'] = 'Multiple Match'
                result_row['match_index'] = match_idx
                result_row['match_count'] = len(matches)
                
                # 标记最佳匹配（MW差值最小的）
                if highlight_best_match:
                    result_row['is_best_match'] = (match_idx == 1)
                
                results.append(result_row)
        
        # 进度显示
        if (idx + 1) % 50 == 0:
            print(f"已处理: {idx + 1}/{len(annotation_df)} 行")
    
    # 创建结果DataFrame
    result_df = pd.DataFrame(results)
    
    # 列排序
    priority_cols = ['original_row_number', 'match_status', 'match_index', 
                    'match_count', 'is_best_match', 'MW_difference', 'abs_MW_difference']
    original_cols = annotation_df.columns.tolist()
    matched_cols = [col for col in result_df.columns if col.startswith('matched_')]
    
    final_cols = []
    for col in priority_cols:
        if col in result_df.columns:
            final_cols.append(col)
    
    final_cols.extend(original_cols)
    final_cols.extend(matched_cols)
    
    # 移除重复
    final_cols = list(dict.fromkeys(final_cols))
    result_df = result_df[final_cols]
    
    # 保存
    print("\n正在保存结果...")
    try:
        # 使用ExcelWriter来应用一些格式
        with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
            result_df.to_excel(writer, index=False, sheet_name='Matched Results')
            
            # 获取工作表
            workbook = writer.book
            worksheet = writer.sheets['Matched Results']
            
            # 调整列宽
            for column in worksheet.columns:
                max_length = 0
                column = [cell for cell in column]
                for cell in column:
                    try:
                        if len(str(cell.value)) > max_length:
                            max_length = len(cell.value)
                    except:
                        pass
                adjusted_width = min(max_length + 2, 50)
                worksheet.column_dimensions[column[0].column_letter].width = adjusted_width
        
        print(f"✓ 结果已保存到: {output_file}")
    except Exception as e:
        print(f"错误: 保存文件失败 - {e}")
        return
    
    # 统计信息
    print("\n" + "=" * 60)
    print("匹配统计结果")
    print("=" * 60)
    print(f"原始待注释行数:     {stats['total_rows']}")
    print(f"输出结果行数:       {len(result_df)}")
    print(f"\n匹配情况分布:")
    print(f"  - 无匹配:         {stats['no_match']} 行 ({stats['no_match']/stats['total_rows']*100:.1f}%)")
    print(f"  - 单个匹配:       {stats['single_match']} 行 ({stats['single_match']/stats['total_rows']*100:.1f}%)")
    print(f"  - 多个匹配:       {stats['multiple_match']} 行 ({stats['multiple_match']/stats['total_rows']*100:.1f}%)")
    
    if stats['multiple_match'] > 0:
        multi_df = result_df[result_df['match_status'] == 'Multiple Match']
        print(f"\n多匹配统计:")
        print(f"  总共扩展出 {len(multi_df)} 行")
        if highlight_best_match:
            best_matches = sum(multi_df['is_best_match'] == True)
            print(f"  其中最佳匹配: {best_matches} 行")
    
    print("=" * 60)
    print(f"处理完成! 时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)
    
    return result_df


# 使用增强版
if __name__ == "__main__":
    result = annotate_all_matches_advanced(
        raw_file='raw.xlsx',
        annotation_file='for_annotation.xlsx',
        output_file='annotated_result_enhanced.xlsx',
        mw_tolerance=0.03,
        add_original_index=True,      # 添加原始行号
        highlight_best_match=True     # 标记最佳匹配
    )
