from Scripts.tools import *


def main():
    smiles = input("Введите SMILES: ")
    generator = MoleculeConformerGenerator(smiles)
    generator.generate_and_save()


if __name__ == "__main__":
    main()