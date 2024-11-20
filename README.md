# Computational-Chemistry

This project converts a SMILES string into a Gaussian input file (.gjf) for geometry optimization and frequency calculations.
It uses RDKit for molecule generation and a custom Gaussian input generator.

## Installation

1. **Prerequisites:**
   * Docker: Ensure Docker Desktop is installed and running on your system.
   * Git: Install Git.
   * Python 3.10 or higher (for local testing, optional)


2. **Clone the Repository:**

   ```bash
   git clone https://github.com/BorisGrechkin/Computational-Chemistry.git
   cd Computational-Chemistry
   ```

3. **Build the Docker Image:**

   ```
   docker build -t smiles_to-gjf:1.0 .
   ```

4. **Usage**
   
   ```bash
   docker run --rm -v <path_to_results>:/Results smiles_to-gjf:1.0 --input_smiles "<smiles_string>"
   ```

*Replace <path_to_results> with your actual absolute path. Always use absolute paths for the volume mount.

Output Files: The output will be generated inside the /Results directory on your host machine.
.mol files will be generated in /Coordinates 
.gjf files will be generated in /GJF_files

# XYZ_to_GJF Converter

Скрипт для конвертации координат атомов в молекуле из формата `.mol` (или `.xyz`) в формат `.gjf`. Этот формат необходим для работы с программой **Gaussian 09W**.

---

## Описание

Данный скрипт используется для преобразования файлов с координатами атомов в молекуле из формата `.mol` (он же `.xyz`) в `.gjf`. Это необходимо для дальнейших расчетов в программе **Gaussian 09W**. 

Скрипт позволяет задавать параметры для расчета: метод, заряд молекулы, мультиплетность, количество процессоров и другие настройки.

---

## Установка

Инструкция по использованию
Шаг 1: Подготовка координат
Откройте папку Coordinates. В ней находятся файлы с расширением .mol, содержащие координаты атомов молекулы. Эти файлы понадобятся для указания пути при запуске скрипта.

Примечание: Расширения .mol и .xyz идентичны в данном случае.

Шаг 2: Запуск скрипта
Откройте скрипт XYZ_to_gjf.py. Для его запуска укажите параметры, определяющие цель расчета. Например:

Метод: B3LYP
Базисный набор (если требуется): укажите отдельно.
Заряд молекулы: 0
Мультиплетность: 1
Количество процессоров: 8
Запуск выполняется командой:

'''bash
python Get_gjf\\XYZ_to_gjf.py --xyz_files "Coordinates\\C1=CC=CC=C1_1.mol" --method B3LYP --charge 0 --multiplicity 1 --nproc 8
'''
Шаг 3: Результаты
После выполнения скрипта файлы с расширением .gjf и указанными параметрами будут сохранены в папке Results. Теперь эти файлы можно использовать в программе Gaussian 09W для расчетов.
