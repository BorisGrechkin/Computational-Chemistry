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
   docker build -t smiles_to-gif:1.0 .
   ```

4. **Usage**
   
   ```bash
   docker run --rm -v <path_to_results>:/Results smiles_to-gif:1.0 --input_smiles "<smiles_string>"
   ```

*Replace /path/to/my/results with your actual absolute path. Always use absolute paths for the volume mount.

Output Files: The output will be generated inside the /Results directory on your host machine.
.mol files will be generated in /Coordinates 
.gjf files will be generated in /GIF_files
