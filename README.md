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

A script for converting atomic coordinates in a molecule from the `.mol` (or `.xyz`) format to the `.gjf` format. This format is essential for performing calculations in the **Gaussian 09W** program.

---

## Description

This script is used to transform files containing atomic coordinates in a molecule from the `.mol` format (which is equivalent to `.xyz`) into `.gjf`. This conversion is necessary for further computations in **Gaussian 09W**.

The script allows you to specify calculation parameters, including the method, molecular charge, multiplicity, the number of processors, and more.

---

## Installation

---

## Instructions for Use

**Step 1: Preparing the coordinates**  
Open the **Coordinates** folder. It contains `.mol` files with atomic coordinates of molecules. These files are needed to specify the input path when running the script.

Note: The `.mol` and `.xyz` extensions are interchangeable in this context.

---

**Step 2: Running the script**  
Open the `XYZ_to_gjf.py` script. Specify the parameters that determine the calculation's goal. For example:

* Method: `B3LYP`
* Basis set (if required): specify separately.
* Molecular charge: `0`
* Multiplicity: `1`
* Number of processors: `8`

The script can be run using the following command:

```bash
python Scripts\\XYZ_to_gjf.py --xyz_files "Results\\Coordinates\\C1=CC=CC=C1_1.mol" --level_of_theory B3LYP/6-31G**  --method "Opt Freq"  --charge 0 --multiplicity 1 --nproc 4 
```

Step 3: Results

After running the script, .gjf files with the specified parameters will be saved in the Results folder. These files are now ready to be used in the Gaussian 09W program for further calculations
