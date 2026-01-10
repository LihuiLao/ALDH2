import pandas as pd
import os
from datetime import datetime

def annotate_all_matches_advanced(raw_file, annotation_file, output_file, 
                                  mw_tolerance=0.03, 
                                  add_original_index=True,
                                  highlight_best_match=True):
    """
    ：
    
    :
        add_original_index: 
        highlight_best_match: （MW）
    """
    
    print("=" * 60)
    print("Excel  ()")
    print("=" * 60)
    
    # 
    if not os.path.exists(raw_file):
        print(f":  {raw_file}")
        return
    if not os.path.exists(annotation_file):
        print(f":  {annotation_file}")
        return
    
    # 
    print(f"\n...")
    try:
        raw_df = pd.read_excel(raw_file)
        print(f"✓  {raw_file}: {len(raw_df)} , {len(raw_df.columns)} ")
        
        annotation_df = pd.read_excel(annotation_file)
        print(f"✓  {annotation_file}: {len(annotation_df)} , {len(annotation_df.columns)} ")
    except Exception as e:
        print(f":  - {e}")
        return
    
    # MW
    if 'MW' not in raw_df.columns or 'MW' not in annotation_df.columns:
        print(":  'MW' ")
        return
    
    print(f"\n...")
    print("-" * 60)
    
    results = []
    stats = {
        'total_rows': len(annotation_df),
        'no_match': 0,
        'single_match': 0,
        'multiple_match': 0,
    }
    
    # 
    for idx, row in annotation_df.iterrows():
        mw_value = row['MW']
        
        # 
        matches = raw_df[
            (raw_df['MW'] >= mw_value - mw_tolerance) & 
            (raw_df['MW'] <= mw_value + mw_tolerance)
        ].copy()
        
        # MW
        if not matches.empty:
            matches['abs_mw_diff'] = abs(matches['MW'] - mw_value)
            matches = matches.sort_values('abs_mw_diff')
        
        if matches.empty:
            # 
            stats['no_match'] += 1
            result_row = row.to_dict()
            if add_original_index:
                result_row['original_row_number'] = idx + 2  # Excel2（）
            result_row['match_status'] = 'No Match'
            result_row['match_index'] = None
            result_row['match_count'] = 0
            result_row['is_best_match'] = None
            results.append(result_row)
            
        elif len(matches) == 1:
            # 
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
            # 
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
                
                # （MW）
                if highlight_best_match:
                    result_row['is_best_match'] = (match_idx == 1)
                
                results.append(result_row)
        
        # 
        if (idx + 1) % 50 == 0:
            print(f": {idx + 1}/{len(annotation_df)} ")
    
    # DataFrame
    result_df = pd.DataFrame(results)
    
    # 
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
    
    # 
    final_cols = list(dict.fromkeys(final_cols))
    result_df = result_df[final_cols]
    
    # 
    print("\n...")
    try:
        # ExcelWriter
        with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
            result_df.to_excel(writer, index=False, sheet_name='Matched Results')
            
            # 
            workbook = writer.book
            worksheet = writer.sheets['Matched Results']
            
            # 
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
        
        print(f"✓ : {output_file}")
    except Exception as e:
        print(f":  - {e}")
        return
    
    # 
    print("\n" + "=" * 60)
    print("")
    print("=" * 60)
    print(f":     {stats['total_rows']}")
    print(f":       {len(result_df)}")
    print(f"\n:")
    print(f"  - :         {stats['no_match']}  ({stats['no_match']/stats['total_rows']*100:.1f}%)")
    print(f"  - :       {stats['single_match']}  ({stats['single_match']/stats['total_rows']*100:.1f}%)")
    print(f"  - :       {stats['multiple_match']}  ({stats['multiple_match']/stats['total_rows']*100:.1f}%)")
    
    if stats['multiple_match'] > 0:
        multi_df = result_df[result_df['match_status'] == 'Multiple Match']
        print(f"\n:")
        print(f"   {len(multi_df)} ")
        if highlight_best_match:
            best_matches = sum(multi_df['is_best_match'] == True)
            print(f"  : {best_matches} ")
    
    print("=" * 60)
    print(f"! : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)
    
    return result_df


# 
if __name__ == "__main__":
    result = annotate_all_matches_advanced(
        raw_file='raw.xlsx',
        annotation_file='for_annotation.xlsx',
        output_file='annotated_result_enhanced.xlsx',
        mw_tolerance=0.03,
        add_original_index=True,      # 
        highlight_best_match=True     # 
    )
