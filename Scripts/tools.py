import os
from rdkit import Chem
from rdkit.Chem import AllChem


class MoleculeConformerGenerator:
    """
    Класс для генерации и оптимизации конформеров молекулы на основе SMILES-строки.

    Атрибуты:
        smiles (str): SMILES-строка молекулы.
        num_conformers (int): Количество генерируемых конформеров. По умолчанию 10.
        output_dir (str): Папка для сохранения конформеров. По умолчанию "Coordinates".
        mol (RDKit Mol): Молекула, полученная из SMILES.
        mol_with_h (RDKit Mol): Молекула с добавленными атомами водорода.
        conformer_ids (list): Список идентификаторов сгенерированных конформеров.
        energies (list): Список энергий конформеров.

    Методы:
        __init__(smiles, num_conformers, output_dir): Инициализирует объект с заданными параметрами.
        get_molecule_from_smiles(): Преобразует SMILES в молекулу.
        add_hydrogens(): Добавляет атомы водорода в молекулу.
        generate_conformers(): Генерирует несколько конформеров молекулы.
        optimize_conformers(): Оптимизирует конформеры и вычисляет их энергии.
        save_conformers(): Сохраняет топ-конформеры в файлы.
        generate_and_save(): Полный процесс генерации и сохранения конформеров.
    """

    def __init__(self, smiles: str, num_conformers: int = 10, output_dir="Results/Coordinates"):
        """
        Инициализирует объект MoleculeConformerGenerator с заданными параметрами.

        Параметры:
            smiles (str): SMILES-строка молекулы.
            num_conformers (int): Количество генерируемых конформеров (по умолчанию 10).
            output_dir (str): Папка для сохранения результатов (по умолчанию "Coordinates").
        """
        self.smiles = smiles
        self.num_conformers = num_conformers
        self.output_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', output_dir))
        self.mol = None
        self.mol_with_h = None
        self.conformer_ids = []
        self.energies = []

    def get_molecule_from_smiles(self):
        """
        Преобразует SMILES-строку в молекулу.

        Возвращает:
            RDKit Mol: Молекула, полученная из SMILES.

        Исключения:
            ValueError: Если SMILES-строка некорректна.
        """
        self.mol = Chem.MolFromSmiles(self.smiles)
        if self.mol is None:
            raise ValueError("Невозможно преобразовать SMILES в молекулу. Проверьте ввод.")
        return self.mol

    # TODO: (thegade) Посмотреть как по близости геометрии фильтровать.
    # можно попросить катю разобраться как это делается, какие надо оставлять при такой фильтрации.
    # добавить метод фильтрации по геометрии

    def add_hydrogens(self):
        """
        Добавляет атомы водорода в молекулу.

        Возвращает:
            RDKit Mol: Молекула с добавленными атомами водорода.
        """
        self.mol_with_h = Chem.AddHs(self.mol)
        return self.mol_with_h

    def generate_conformers(self):
        """
        Генерирует несколько конформеров молекулы.

        Возвращает:
            list: Список идентификаторов сгенерированных конформеров.
        """
        self.conformer_ids = AllChem.EmbedMultipleConfs(self.mol_with_h, numConfs=self.num_conformers, params=AllChem.ETKDG())
        return self.conformer_ids

    def optimize_conformers(self):
        """"
        Оптимизирует конформеры и вычисляет их энергии.

        Возвращает:
            list: Список кортежей с идентификаторами конформеров и их энергиями.
        """
        self.energies = []
        for conf_id in self.conformer_ids:
            AllChem.UFFOptimizeMolecule(self.mol_with_h, confId=conf_id)
            ff = AllChem.UFFGetMoleculeForceField(self.mol_with_h, confId=conf_id)
            energy = ff.CalcEnergy()
            self.energies.append((conf_id, energy))
        return self.energies

    def save_conformers(self):
        """
        Сохраняет топ-3 стабильных конформера в формате .mol, содержащем только координаты атомов.
        Формат сохраняемых данных:
        - атомный символ
        - координаты (x, y, z) с точностью до 4 знаков после запятой
        Печатает:
            Сообщение о том, что топ-3 стабильных конформеров сохранены в указанной папке.
        """
        os.makedirs(self.output_dir, exist_ok=True)
        sorted_energies = sorted(self.energies, key=lambda x: x[1])
        top_3_conformers = sorted_energies[:3]

        for conf_id, energy in top_3_conformers:
            file_path = os.path.join(self.output_dir, f"{self.smiles}_{conf_id}.mol")
            conf = self.mol_with_h.GetConformer(conf_id)

            with open(file_path, "w") as file:
                for atom in self.mol_with_h.GetAtoms():
                    atom_symbol = atom.GetSymbol()
                    pos = conf.GetAtomPosition(atom.GetIdx())
                    file.write(f"{atom_symbol:>1}  {pos.x:7.4f}  {pos.y:7.4f}  {pos.z:7.4f}\n")
            print(f"Generated: {self.smiles}_{conf_id}.mol")

    def generate_and_save(self):
        """
        Полный процесс генерации, оптимизации и сохранения конформеров.

        Печатает:
            Сообщение о том, что топ-3 стабильных конформеров сохранены в указанной папке.
        """

        # Тут нужно будет добавить флаг, по какому методу фильтровать
        # -по энергии
        # -по близости геометрии

        self.get_molecule_from_smiles()
        self.add_hydrogens()
        self.generate_conformers()
        self.optimize_conformers()
        self.save_conformers()

class GaussianInputGenerator:
    """
    Generates Gaussian input files (.gjf) from XYZ coordinate files.

    This class provides a convenient way to create Gaussian input files with customizable parameters,
    making it easier to run various quantum chemical calculations.  It handles reading coordinates from
    XYZ files and writing the formatted input to .gjf files.
    """

    def __init__(self, xyz_filepath = 'Results/Coordinates', output_dir = 'Results/GJF_files',
                 method='Opt Freq SCRF=(CPCM,Solvent=Water)', comment='Here could be your advertisement',
                 charge=0, multiplicity=1, title='Molecule', level_of_theory='PM3', nproc=4):
        """
        Initializes GaussianInputGenerator with default parameters.

        Args:
            xyz_filepath (str): Path to the input file with XYZ files.
            output_dir (str): Directory to save the generated .gjf file.
            method (str): Gaussian calculation method (default: 'Opt Freq SCRF=(CPCM,Solvent=Water)').
            comment (str): Comment to be included in the Gaussian input file (default: 'Here could be your advertisement').
            charge (int): Charge of the molecule (default: 0).
            multiplicity (int): Spin multiplicity of the molecule (default: 1).
            title (str): Title for the Gaussian input file (default: 'Molecule').
            level_of_theory (str): Level of theory specification (default: 'PM3').
            nproc (int): Number of processors to be used in the calculation (default: 4).
        """
        self.xyz_filepath = xyz_filepath
        self.output_dir = output_dir
        self.method = method
        self.comment = comment
        self.charge = charge
        self.multiplicity = multiplicity
        self.title = title
        self.level_of_theory = level_of_theory
        self.nproc = nproc

        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

    def generate_gjf(self, coords, filename):
        """
        Generates a Gaussian input file (.gjf).

        Args:
            coords (str): String containing the XYZ coordinates of the molecule.
            filename (str): Path to the output .gjf file.
        """
        with open(filename, 'w') as f:
            f.write(f'%NProcShared={self.nproc}\n')
            f.write(f'# {self.method} {self.comment} {self.level_of_theory}\n\n')
            f.write(f'{self.title}\n\n')
            f.write(f'{self.charge} {self.multiplicity}\n')
            f.write(coords)
            f.write('\n')

    def read_coords_from_xyz(self, file_path):
        """
        Reads coordinates from an XYZ file.

        Args:
            file_path (str): Path to the XYZ file.

        Returns:
            str: String containing the XYZ coordinates, or None if the file is not found.
        """
        try:
            with open(file_path, 'r') as file:
                coords = file.read()
            return coords
        except FileNotFoundError:
            print(f"Error: File not found: {file_path}")
            return None

    def process_xyz_file(self):
        """
        Processes a single XYZ file and generates a Gaussian input file.
        """
        for file in os.listdir(self.xyz_filepath):

            if file[-4:] != '.mol':
                pass
            else:
                coords = self.read_coords_from_xyz(self.xyz_filepath + '/' + file)
                if coords is not None:
                    filename = os.path.join(self.output_dir, os.path.splitext(os.path.basename(file))[0] + ".gjf")
                    self.generate_gjf(coords, filename)
                    print(f'Generated {filename} from file: {file}')
