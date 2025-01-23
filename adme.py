import os
import csv
import shutil
from rdkit import Chem
from rdkit.Chem import Crippen, Lipinski, Descriptors

input_dir = "split_sdf_files"
output_csv = "adme.csv"
adme_pass_dir = "adme-pass"
adme_failed_dir = "adme-failed"
summary_file = "summary-adme.txt"

# Create output directories if they don't exist
os.makedirs(adme_pass_dir, exist_ok=True)
os.makedirs(adme_failed_dir, exist_ok=True)

# Initialize counters and lists for summary
pass_count = 0
fail_count = 0
violate_0_criteria = 0
violate_1_criteria = 0
violate_2_criteria = 0
violate_3_criteria = 0
violate_4_criteria = 0

# Function to check Lipinski's Rule of Five
def check_lipinski(mol):
    mw = Descriptors.MolWt(mol)  # Correct function for molecular weight
    logP = Crippen.MolLogP(mol)
    hdonors = Lipinski.NumHDonors(mol)
    hacceptors = Lipinski.NumHAcceptors(mol)

    failed_criteria = []
    if mw > 500:
        failed_criteria.append("Molecular Weight > 500")
    if logP > 5:
        failed_criteria.append("logP > 5")
    if hdonors > 5:
        failed_criteria.append("NumHDonors > 5")
    if hacceptors > 10:
        failed_criteria.append("NumHAcceptors > 10")

    result = {
        'MolecularWeight': mw,
        'logP': logP,
        'NumHDonors': hdonors,
        'NumHAcceptors': hacceptors,
        'LipinskiPassed': len(failed_criteria) < 2,  # Pass if fewer than 2 criteria fail
        'FailedCriteria': ", ".join(failed_criteria) if failed_criteria else "None",
        'NumViolations': len(failed_criteria)
    }
    return result

# Create CSV file and write header
with open(output_csv, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['FileName', 'MoleculeName', 'MolecularWeight', 'logP', 'NumHDonors', 'NumHAcceptors', 
                     'LipinskiPassed', 'FailedCriteria', 'NumViolations'])

    # Iterate over SDF files
    for filename in os.listdir(input_dir):
        if filename.endswith(".sdf"):
            input_path = os.path.join(input_dir, filename)

            # Read molecules from SDF
            suppl = Chem.SDMolSupplier(input_path)

            for i, mol in enumerate(suppl):
                if mol is not None:
                    result = check_lipinski(mol)
                    molecule_name = mol.GetProp('_Name') if mol.HasProp('_Name') else f'Molecule_{i+1}'

                    # Write data to CSV
                    writer.writerow([filename, molecule_name, result['MolecularWeight'], result['logP'], 
                                     result['NumHDonors'], result['NumHAcceptors'], 
                                     'Passed' if result['LipinskiPassed'] else 'Failed',
                                     result['FailedCriteria'], result['NumViolations']])

                    # Track statistics for summary
                    num_violations = result['NumViolations']
                    if num_violations == 0:
                        violate_0_criteria += 1
                    elif num_violations == 1:
                        violate_1_criteria += 1

                    if not result['LipinskiPassed']:  # If the compound failed Lipinski's rules
                        fail_count += 1
                        if num_violations == 2:
                            violate_2_criteria += 1
                        elif num_violations == 3:
                            violate_3_criteria += 1
                        elif num_violations == 4:
                            violate_4_criteria += 1
                        shutil.copy(input_path, os.path.join(adme_failed_dir, filename))
                    else:
                        pass_count += 1
                        shutil.copy(input_path, os.path.join(adme_pass_dir, filename))

# Write the summary to a text file
with open(summary_file, mode='w') as summary:
    summary.write(f"Total Compounds Processed: {pass_count + fail_count}\n")
    summary.write(f"Compounds that Passed Lipinski's Rules: {pass_count}\n")
    summary.write(f"Compounds that Failed Lipinski's Rules: {fail_count}\n")
    summary.write(f"  - Compounds that Violate 0 Criterion: {violate_0_criteria}\n")
    summary.write(f"  - Compounds that Violate 1 Criterion: {violate_1_criteria}\n")
    summary.write(f"  - Compounds that Violate 2 Criteria: {violate_2_criteria}\n")
    summary.write(f"  - Compounds that Violate 3 Criteria: {violate_3_criteria}\n")
    summary.write(f"  - Compounds that Violate 4 Criteria: {violate_4_criteria}\n")

print(f"Results saved in {output_csv}")
print(f"Compounds that passed Lipinski's rules are in '{adme_pass_dir}' and those that failed are in '{adme_failed_dir}'.")
print(f"Summary saved in {summary_file}")

