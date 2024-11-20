import argparse

from Scripts.tools import *


def main():

    parser = argparse.ArgumentParser(description='Generate .gjf files from SMILES')
    parser.add_argument('--input_smiles', type=str,
                        required=True, help='Input SMILES formula')

    args = parser.parse_args()

    conformer_generator = MoleculeConformerGenerator(args.input_smiles)
    conformer_generator.generate_and_save()
    gif_generator = GaussianInputGenerator()
    gif_generator.process_xyz_file()


if __name__ == "__main__":
    main()