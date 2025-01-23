import os
import argparse

def split_sdf(input_file, output_dir):
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    with open(input_file, 'r') as sdf_file:
        content = sdf_file.read()

    # Split the file by molecule entry (identified by "$$$$")
    molecules = content.split("$$$$\n")
    
    for molecule in molecules:
        if molecule.strip():  # Ignore empty splits
            # Extract the ID from the <npaid> tag
            lines = molecule.splitlines()
            for i, line in enumerate(lines):
                if line.strip().startswith(">  <npaid>"):
                    npaid = lines[i + 1].strip()
                    break
            else:
                print("Warning: No <npaid> field found in one molecule entry.")
                continue
            
            # Write the molecule to a new SDF file
            output_file = os.path.join(output_dir, f"{npaid}.sdf")
            with open(output_file, 'w') as out_file:
                out_file.write(molecule + "$$$$\n")
    
    print(f"Splitting complete. Molecules saved in {output_dir}.")

if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Split an SDF file into individual files based on <npaid>.")
    parser.add_argument("-i", "--input", required=True, help="Input SDF file path")
    parser.add_argument("-o", "--output", default="split_sdf_files", help="Output directory (default: 'split_sdf_files')")
    
    args = parser.parse_args()
    
    # Run the split function
    split_sdf(args.input, args.output)

