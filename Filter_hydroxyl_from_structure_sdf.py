from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski  #  Lipinski 
import pandas as pd

def extract_hydroxyl_compounds(sdf_file='structures.sdf'):
    """
    HMDBSDF
    """
    print("SDF...")
    suppl = Chem.SDMolSupplier(sdf_file)
    
    # 
    hydroxyl_pattern = Chem.MolFromSmarts('[OX2H]')
    
    results = []
    total = 0
    found = 0
    
    for mol in suppl:
        total += 1
        
        if mol is None:
            continue
        
        # 
        if mol.HasSubstructMatch(hydroxyl_pattern):
            found += 1
            
            # 
            num_oh = len(mol.GetSubstructMatches(hydroxyl_pattern))
            
            # 
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
                'HBD': Lipinski.NumHDonors(mol),  # 
                'HBA': Lipinski.NumHAcceptors(mol),  # 
            })
            
            # 
            if found % 500 == 0:
                print(f":  {total} ， {found} ")
    
    print(f"\n！ {total} ， {found} ")
    
    # 
    df = pd.DataFrame(results)
    df.to_csv('hydroxyl_compounds.csv', index=False, encoding='utf-8-sig')
    df.to_excel('hydroxyl_compounds.xlsx', index=False)
    print(f"\n：")
    print(f"  - hydroxyl_compounds.csv")
    print(f"  - hydroxyl_compounds.xlsx")
    
    # 
    save_structures(sdf_file, 'hydroxyl_compounds.sdf')
    
    return df

def save_structures(input_sdf, output_sdf):
    """
    SDF
    """
    print("\n...")
    suppl = Chem.SDMolSupplier(input_sdf)
    writer = Chem.SDWriter(output_sdf)
    hydroxyl_pattern = Chem.MolFromSmarts('[OX2H]')
    
    count = 0
    for mol in suppl:
        if mol is not None and mol.HasSubstructMatch(hydroxyl_pattern):
            writer.write(mol)
            count += 1
    
    writer.close()
    print(f"  - {output_sdf} ( {count} )")

if __name__ == "__main__":
    # 
    df = extract_hydroxyl_compounds('structures.sdf')
    
    # 
    print("\n===  ===")
    print(f": {len(df)}")
    print(f": {df['MW'].mean():.2f}")
    print(f": {df['MW'].min():.2f} - {df['MW'].max():.2f}")
    print(f": {df['Num_OH'].mean():.2f}")
    print(f": {df['Num_OH'].max()}")
    
    # 10
    print("\n=== 10 ===")
    print(df[['HMDB_ID', 'Name', 'Formula', 'Num_OH']].head(10).to_string())
