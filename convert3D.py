import os
from rdkit import Chem
from rdkit.Chem import AllChem
from concurrent.futures import ProcessPoolExecutor, as_completed

# Input and output directories
input_dir = "adme-pass"
output_dir = "out3D"
error_log = "error.txt"
num_cores = 4  # Default to 4 cores

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Function to process each molecule and generate 3D coordinates
def process_molecule(filename, input_path, output_path):
    error_message = None
    try:
        suppl = Chem.SDMolSupplier(input_path)
        writer = Chem.SDWriter(output_path)

        for mol in suppl:
            if mol is None:
                error_message = f"Warning: Could not read molecule in {filename}"
                continue

            # Add hydrogens to the molecule
            mol = Chem.AddHs(mol)

            # Generate 3D coordinates
            AllChem.EmbedMolecule(mol, AllChem.ETKDG())
            AllChem.UFFOptimizeMolecule(mol)
            writer.write(mol)

        writer.close()

    except Exception as e:
        error_message = f"Error processing {filename}: {e}"

    return error_message

# Function to manage multi-core processing
def process_files():
    with open(error_log, "w") as error_file:
        # Create a list to hold all the future tasks
        futures = []
        with ProcessPoolExecutor(max_workers=num_cores) as executor:
            for filename in os.listdir(input_dir):
                if filename.endswith(".sdf"):
                    input_path = os.path.join(input_dir, filename)
                    output_path = os.path.join(output_dir, filename)

                    # Submit tasks to the executor for parallel processing
                    futures.append(executor.submit(process_molecule, filename, input_path, output_path))

            # Process completed tasks and handle any errors
            for future in as_completed(futures):
                error_message = future.result()
                if error_message:
                    error_file.write(error_message + "\n")
                    print(error_message)

    print("3D conversion completed. Files saved in 'out3D'. Check 'error.txt' for errors.")

# Run the function
process_files()

