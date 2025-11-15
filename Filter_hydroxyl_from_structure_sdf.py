from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski  # 添加 Lipinski 导入
import pandas as pd

def extract_hydroxyl_compounds(sdf_file='structures.sdf'):
    """
    从HMDB的SDF文件中提取含羟基的化合物
    """
    print("开始读取SDF文件...")
    suppl = Chem.SDMolSupplier(sdf_file)
    
    # 定义羟基模式
    hydroxyl_pattern = Chem.MolFromSmarts('[OX2H]')
    
    results = []
    total = 0
    found = 0
    
    for mol in suppl:
        total += 1
        
        if mol is None:
            continue
        
        # 检查是否含有羟基
        if mol.HasSubstructMatch(hydroxyl_pattern):
            found += 1
            
            # 计算羟基数量
            num_oh = len(mol.GetSubstructMatches(hydroxyl_pattern))
            
            # 提取属性
            try:
                hmdb_id = mol.GetProp('DATABASE_ID') if mol.HasProp('DATABASE_ID') else 'N/A'
            except:
                hmdb_id = 'N/A'
                
            try:
                name = mol.GetProp('GENERIC_NAME') if mol.HasProp('GENERIC_NAME') else mol.GetProp('_Name')
            except:
                name = 'Unknown'
            
            results.append({
                'HMDB_ID': hmdb_id,
                'Name': name,
                'SMILES': Chem.MolToSmiles(mol),
                'Formula': Chem.rdMolDescriptors.CalcMolFormula(mol),
                'MW': round(Descriptors.MolWt(mol), 2),
                'Num_OH': num_oh,
                'LogP': round(Descriptors.MolLogP(mol), 2),
                'HBD': Lipinski.NumHDonors(mol),  # 氢键供体
                'HBA': Lipinski.NumHAcceptors(mol),  # 氢键受体
            })
            
            # 进度提示
            if found % 500 == 0:
                print(f"进度: 已处理 {total} 个，找到 {found} 个含羟基化合物")
    
    print(f"\n完成！总共处理 {total} 个化合物，找到 {found} 个含羟基化合物")
    
    # 保存结果
    df = pd.DataFrame(results)
    df.to_csv('hydroxyl_compounds.csv', index=False, encoding='utf-8-sig')
    df.to_excel('hydroxyl_compounds.xlsx', index=False)
    print(f"\n结果已保存：")
    print(f"  - hydroxyl_compounds.csv")
    print(f"  - hydroxyl_compounds.xlsx")
    
    # 保存结构文件
    save_structures(sdf_file, 'hydroxyl_compounds.sdf')
    
    return df

def save_structures(input_sdf, output_sdf):
    """
    保存含羟基化合物的结构到新的SDF文件
    """
    print("\n正在保存结构文件...")
    suppl = Chem.SDMolSupplier(input_sdf)
    writer = Chem.SDWriter(output_sdf)
    hydroxyl_pattern = Chem.MolFromSmarts('[OX2H]')
    
    count = 0
    for mol in suppl:
        if mol is not None and mol.HasSubstructMatch(hydroxyl_pattern):
            writer.write(mol)
            count += 1
    
    writer.close()
    print(f"  - {output_sdf} (包含 {count} 个结构)")

if __name__ == "__main__":
    # 直接运行
    df = extract_hydroxyl_compounds('structures.sdf')
    
    # 显示基本统计
    print("\n=== 统计摘要 ===")
    print(f"化合物总数: {len(df)}")
    print(f"平均分子量: {df['MW'].mean():.2f}")
    print(f"分子量范围: {df['MW'].min():.2f} - {df['MW'].max():.2f}")
    print(f"平均羟基数: {df['Num_OH'].mean():.2f}")
    print(f"最多羟基数: {df['Num_OH'].max()}")
    
    # 显示前10个结果
    print("\n=== 前10个化合物 ===")
    print(df[['HMDB_ID', 'Name', 'Formula', 'Num_OH']].head(10).to_string())
