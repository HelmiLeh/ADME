import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem


def detect_delimiter(file_path):
    ext = os.path.splitext(file_path)[1].lower()
    if ext in ['.tsv']:
        return '\t'
    elif ext in ['.csv']:
        return ','
    elif ext in ['.smi']:
        return 'whitespace'
    else:
        # Try autodetect
        with open(file_path, 'r') as f:
            line = f.readline()
            if '\t' in line:
                return '\t'
            elif ',' in line:
                return ','
            else:
                return 'whitespace'


def read_smiles_file(file_path):
    delimiter = detect_delimiter(file_path)

    if delimiter == 'whitespace':
        # Expect: SMILES [optional name/fields]
        df = pd.read_csv(file_path, sep='\s+', header=None)
        df.columns = ['smiles'] + [f'col_{i}' for i in range(1, df.shape[1])]
    else:
        df = pd.read_csv(file_path, sep=delimiter)
        if 'smiles' not in df.columns:
            raise ValueError("No 'smiles' column found in input file.")
    return df


def convert_cxsmiles_to_sdf(input_file, output_file, fail_smi_file='fail.smi', log_file='fail_log.txt'):
    df = read_smiles_file(input_file)

    sdf_writer = Chem.SDWriter(output_file)
    fail_writer = open(fail_smi_file, 'w')
    log_writer = open(log_file, 'w')

    for index, row in df.iterrows():
        smiles = row['smiles']
        mol = Chem.MolFromSmiles(smiles)

        if mol is None:
            fail_writer.write(smiles + '\n')
            log_writer.write(f"Row {index}: Failed to parse SMILES\n")
            continue

        mol = Chem.AddHs(mol)
        result = AllChem.EmbedMolecule(mol, AllChem.ETKDG())

        if result != 0:
            fail_writer.write(smiles + '\n')
            log_writer.write(f"Row {index}: 3D embedding failed\n")
            continue

        try:
            AllChem.UFFOptimizeMolecule(mol)
        except ValueError as e:
            fail_writer.write(smiles + '\n')
            log_writer.write(f"Row {index}: UFF optimization error: {e}\n")
            continue

        # Add metadata
        for col in df.columns:
            if col != 'smiles':
                mol.SetProp(str(col), str(row[col]))

        sdf_writer.write(mol)

    sdf_writer.close()
    fail_writer.close()
    log_writer.close()
    print(f"✅ Conversion completed: {output_file}")
    print(f"⚠️ Failed SMILES saved to: {fail_smi_file}")
    print(f"📄 Log written to: {log_file}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Convert SMILES (.smi/.csv/.tsv) to 3D SDF with optional metadata.")
    parser.add_argument("-i", "--input", required=True, help="Input file (.smi, .csv, .tsv)")
    parser.add_argument("-o", "--output", required=True, help="Output SDF file (.sdf)")
    parser.add_argument("--fail", default="fail.smi", help="Output file for failed SMILES [default: fail.smi]")
    parser.add_argument("--log", default="fail_log.txt", help="Log file for errors [default: fail_log.txt]")

    args = parser.parse_args()
    convert_cxsmiles_to_sdf(args.input, args.output, args.fail, args.log)

