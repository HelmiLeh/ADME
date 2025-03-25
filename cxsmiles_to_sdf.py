import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

def convert_cxsmiles_to_sdf(input_file, output_file):
    # Read input CXSMILES file (tab-separated)
    df = pd.read_csv(input_file, sep='\t')
    
    # Check if the first column contains SMILES strings
    if 'smiles' not in df.columns:
        raise ValueError("Error: First column must contain SMILES strings with header 'smiles'.")

    # Create an SDF writer
    writer = Chem.SDWriter(output_file)

    for index, row in df.iterrows():
        smiles = row['smiles']
        mol = Chem.MolFromSmiles(smiles)  # Convert SMILES to RDKit molecule
        
        if mol:
            mol = Chem.AddHs(mol)  # Add hydrogens
            AllChem.EmbedMolecule(mol, AllChem.ETKDG())  # Generate 3D coordinates
            AllChem.UFFOptimizeMolecule(mol)  # Optimize structure
            
            # Write molecule to SDF file
            writer.write(mol)
    
    writer.close()
    print(f"Conversion completed! Saved as {output_file}")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Convert CXSMILES to SDF")
    parser.add_argument("-i", "--input", required=True, help="Input CXSMILES file (.cxsmiles)")
    parser.add_argument("-o", "--output", required=True, help="Output SDF file (.sdf)")

    args = parser.parse_args()
    convert_cxsmiles_to_sdf(args.input, args.output)

