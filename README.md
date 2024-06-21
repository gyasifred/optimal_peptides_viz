# optimal_peptides_viz

## Overview
The `optimal_peptides_viz` repository contains a Python package for generating boxplot visualizations of optimal samples based on user-provided datasets. This README.md file serves as a guide for users on how to use the package and provides information about the repository structure.

# Dataset Preparation
Before using the `optimal_peptides` package, ensure that your dataset meets the following requirements:

1. **CSV/TSV File Format:** The dataset should be in either CSV (Comma-Separated Values) or TSV (Tab-Separated Values) format.
2. **Header Naming Convention:** The header of your dataset should follow this naming convention: `sample-sampletype_library_capture protein`

For Example:
```bash 
2022-CW2_PM1_pG
```
The first column must contain the Sequence name.

## Repository Structure
The repository includes the following files:

- `README.md`: This file contains an overview of the `optimal_peptides_viz` package and instructions for users.
- `optimal.py`: The Python script that implements the functionality for generating boxplot visualizations of optimal samples.
- `requirements.txt`: A text file listing the required Python packages and their versions for running the `optimal_peptides_viz` package.

## Usage
To use the `optimal_peptides_viz` package, follow these steps:

1. **Clone the repository to your local machine**
   ```bash
   git clone https://github.com/gyasifred/optimal_peptides_viz.git
   ```
2. **Navigate to the repository directory**
   ```bash
   cd optimal_peptides_viz
   ```
3. **Create a virtual environment (optional but recommended)**
   ```bash
   python -m venv venv
   ```
4.  **Activate the Virtual Environment:**
  - On Windows:
    ```
    .\venv\Scripts\activate
    ```
  - On macOS/Linux:
    ```
    source venv/bin/activate
    ```
5. **Install the required packages**
     ```bash
     pip install -r requirements.txt
     
6. Run the `optimal.py`script with your desired parameters. Use the `--help` option to see available command-line arguments
   ```bash
   python optimal.py --help
   ```

7. **Example Command:**
```bash
python optimal.py my_dataset.csv my_plot --sampletypeinitial M MW MS --peptide_threshold 10 --figsize 8 6 --verbose --operator diff --reverse --log2
```

This command will generate a boxplot visualization of enriched samples based on the provided dataset with specified parameters.

## Additional Notes
- For any issues or questions, please refer to the package documentation or contact the package maintainers  **gyasifred@gmail.com** or contact directly `Emmanuel Agyei Frimpong`**kojomiguel@gmail.com** 


  
