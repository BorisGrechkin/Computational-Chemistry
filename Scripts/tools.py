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

    def __init__(self, smiles: str, num_conformers: int = 10, output_dir="Coordinates"):
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
            print(f"Конформер {conf_id} сохранен в формате .mol: {file_path}")

    def generate_and_save(self):
        """
        Полный процесс генерации, оптимизации и сохранения конформеров.

        Печатает:
            Сообщение о том, что топ-3 стабильных конформеров сохранены в указанной папке.
        """
        self.get_molecule_from_smiles()
        self.add_hydrogens()
        self.generate_conformers()
        self.optimize_conformers()
        self.save_conformers()