
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
import argparse
import os

# Function to generate .gjf file
def generate_gjf(coords, filename, method, charge, multiplicity, title, level_of_theory, nproc, comment):
    with open(filename, 'w') as f:
        # Recording Headers
        f.write(f'%NProcShared={nproc}\n')
        f.write(f'# {method} {comment} {level_of_theory}\n\n')
        f.write(f'{title}\n\n')
        f.write(f'{charge} {multiplicity}\n')
        # Recording coordinates
        f.write(coords)
        f.write('\n')

# Main block using argparse
def main():
    parser = argparse.ArgumentParser(description='Generate .gjf files from SMILES strings with customizable parameters.')
    parser.add_argument('--xyz_files', type=str, nargs='+', required=True, help='Enter your path to xyz file')
    parser.add_argument('--method', type=str, default='Opt Freq SCRF=(CPCM,Solvent=Water)', help='Method for calculation (default: Opt Freq).')
    parser.add_argument('--comment', type=str, default='Here could be your advertisement', help='Just comment.')
    parser.add_argument('--charge', type=int, default=0, help='Charge of the molecule (default: 0).')
    parser.add_argument('--multiplicity', type=int, default=1, help='Multiplicity of the molecule (default: 1).')
    parser.add_argument('--title', type=str, default='Molecule', help='Title for the .gjf file (default: Molecule).')
    parser.add_argument('--level_of_theory', type=str, default='PM3', help='Additional parameters for calculation (default: PM3)).')
    parser.add_argument('--nproc', type=int, default=4, help='Number of processors to use (default: 4).')
    

    args = parser.parse_args()

    def read_coords_from_xyz(file_path):
        with open(file_path, 'r') as file:
            coords = file.read()
        return coords

# Generate .gjf files
    for i, xyz_file in enumerate(args.xyz_files):
        # Check if the file exists
        if not os.path.exists(xyz_file):
            print(f"File not found: {xyz_file}")
            continue

    # Reading coordinates from .xyz file
    coords = read_coords_from_xyz(xyz_file)

    # Generate .gjf file name
    filename = f"molecule_{i+1}.gjf"
    title = f'{args.title} {i+1}'  # Adding an index to the title

    # Calling the generate_gjf function with coordinates from a file
    generate_gjf(
        coords=coords,
        filename=filename,
        method=args.method,
        comment=args.comment,
        charge=args.charge,
        multiplicity=args.multiplicity,
        title=title,
        level_of_theory=args.level_of_theory,
        nproc=args.nproc,
    )

    print(f'Generated {filename} from file: {xyz_file}')
    
if __name__ == '__main__':
    main()    



