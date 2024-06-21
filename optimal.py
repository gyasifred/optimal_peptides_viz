import argparse  # Importing the argparse module for command-line argument parsing
import numpy as np
import pandas as pd  # Importing pandas for data manipulation
import matplotlib.pyplot as plt  # Importing matplotlib for plotting
import seaborn as sns  # Importing seaborn for enhanced visualizations
import os  # Importing os for operating system-related functionalities

# Ignore runtime warnings to keep the output clean
import warnings
warnings.filterwarnings("ignore")

if __name__ == "__main__":
    """ This script generates a boxplot visualization of enriched samples. """

    # Initialize the argument parser
    parser = argparse.ArgumentParser(
        description="Generates a Boxplot Visualization of Enriched Peptides")

    # Add command-line arguments
    parser.add_argument("file", type=str, help="Input TSV or CSV file path")
    parser.add_argument("outputName", type=str, help="Name to save the output file")
    parser.add_argument("output_path", type=str, nargs='?', default=os.getcwd(),
                        help="Output PDF plot file path. Default is the current directory.")
    parser.add_argument("--sampletypeinitial", type=str, nargs='+', default=["M", "MW", "MS"],
                        help="List of initial sample types to include in the plot")
    parser.add_argument("--peptide_threshold", type=float, default=10,
                        help="Threshold for filtering Enriched Peptides")
    parser.add_argument("--figsize", type=int, nargs=2, default=(6, 4),
                        help="Figure size (width and height) for the plot")
    parser.add_argument("--verbose", action="store_true", default=False,
                        help="Enable verbose output. Default is False.")
    parser.add_argument("--operator", type=str, choices=["minus", "diff", "division", "multiplication"],
                        default="minus", help="Operation to perform between sample types.")
    parser.add_argument("--reverse", action="store_true", default=False,
                        help="Reverse the operation order (e.g., col2 - col1 instead of col1 - col2).")
    parser.add_argument("--log2", action="store_true", default=False,
                        help="Apply log base 2 transformation to the data.")

    # Parse the command-line arguments
    args = parser.parse_args()

    if args.verbose:  # Check if verbose mode is enabled
        print("Verbose mode enabled. Processing data...")

    try:
        # Check if the input file exists
        if not os.path.isfile(args.file):
            raise FileNotFoundError(f"File '{args.file}' does not exist.")

        if args.verbose:  # If verbose mode is enabled, print reading status
            print(f"Reading data from {args.file}...")

        # Read data from TSV or CSV file based on the file extension
        if args.file.lower().endswith('.tsv'):
            df = pd.read_csv(args.file, sep='\t', index_col=0, encoding="utf-8")  # Read TSV file
        elif args.file.lower().endswith('.csv'):
            df = pd.read_csv(args.file, index_col=0, encoding="utf-8")  # Read CSV file
        else:
            raise ValueError("Input file must be a TSV or CSV file.")

        # Check column types and convert to float if needed
        for col in df.columns:
            if not np.issubdtype(df[col].dtype, np.number):
                df[col] = df[col].astype(float)

        # Apply log base 10 transformation if --log2 is specified
        if args.log2:
            df = np.log2(df + 1)  # Adding 1 to avoid log(0)
            args.peptide_threshold = np.log2(args.peptide_threshold + 1)

    except Exception as e:
        print(f"Error reading or processing file: {e}")
        exit(1)  # Exit the script if an error occurs

    if args.verbose:  # If verbose mode is enabled, print processing status
        print("Data successfully loaded.")
        print("Processing and analyzing data...")

    columns = list(df.columns)  # Get the list of column names
    samples = {}  # Initialize a dictionary to store sample data

    # Group columns based on their prefixes
    for item in columns:
        prefix = item.split('-')[0]
        samples.setdefault(prefix, []).append(item)

    filtered_columns = []  # Initialize a list to store filtered columns

    # Filter columns based on sample type initials provided
    for sample in samples.values():
        columns = [col for col in sample if col.split("-")[1].split("_")[0] in args.sampletypeinitial]
        filtered_columns.append(columns)

    enriched_peptides_df = {}  # Initialize a dictionary for enriched peptides data

    # Process data to filter enriched peptides based on threshold
    for sample in filtered_columns:
        enriched_peptides = set()
        for sample_type in sample:
            peptides = df.index[df[sample_type] >= args.peptide_threshold]
            enriched_peptides.update(peptides)
        enriched_peptides_df[sample[0].split("-")[0]] = df.loc[list(enriched_peptides), sample]

    sample_type_diff_dfs = {}  # Initialize a dictionary for sample type differences

    # Perform operations between sample types based on the selected operator
    for key, df in enriched_peptides_df.items():
        cols = list(df.columns)
        for i in range(len(cols)):
            for j in range(i + 1, len(cols)):
                col1 = cols[i]
                col2 = cols[j]
                if args.operator == "minus" or args.operator == "diff":
                    # Reverse operation order if --reverse flag is enabled
                    if args.reverse:
                        diff_col_name = f'{col2.split("-")[1].split("_")[0]}_{col1.split("-")[1].split("_")[0]}'
                        df[diff_col_name] = df[col2] - df[col1]
                    else:
                        diff_col_name = f'{col1.split("-")[1].split("_")[0]}_{col2.split("-")[1].split("_")[0]}'
                        df[diff_col_name] = df[col1] - df[col2]
                elif args.operator == "division":
                    # Reverse operation order if --reverse flag is enabled
                    if args.reverse:
                        div_col_name = f'{col2.split("-")[1].split("_")[0]}_{col1.split("-")[1].split("_")[0]}_division'
                        df[div_col_name] = df[col2] / df[col1]
                    else:
                        div_col_name = f'{col1.split("-")[1].split("_")[0]}_{col2.split("-")[1].split("_")[0]}_division'
                        df[div_col_name] = df[col1] / df[col2]
                elif args.operator == "multiplication":
                    # Reverse operation order if --reverse flag is enabled
                    if args.reverse:
                        mul_col_name = f'{col2.split("-")[1].split("_")[0]}_{col1.split("-")[1].split("_")[0]}_multiplication'
                        df[mul_col_name] = df[col2] * df[col1]
                    else:
                        mul_col_name = f'{col1.split("-")[1].split("_")[0]}_{col2.split("-")[1].split("_")[0]}_multiplication'
                        df[mul_col_name] = df[col1] * df[col2]
        # Drop original columns after operations
        df.drop(columns=cols, inplace=True)
        sample_type_diff_dfs[key] = df

    dfs = []  # Initialize a list to store dataframes for plotting

    # Prepare dataframes for plotting
    for sample, df in sample_type_diff_dfs.items():
        index_col_name = df.index.name
        df.reset_index(inplace=True)
        df_melted = pd.melt(df, id_vars=[index_col_name], var_name='Sample Type', value_name='Z-scores')
        df_melted['samples'] = sample
        dfs.append(df_melted)

    # Combine dataframes for plotting
    combined_df = pd.concat(dfs, ignore_index=True)

    # Create a figure and axis
    fig, ax = plt.subplots(figsize=args.figsize, facecolor="white")

    # Plot boxplot with specified parameters
    sns.boxplot(x=combined_df['samples'], y=combined_df['Z-scores'],
                hue=combined_df['Sample Type'], ax=ax, width=0.8)


    ax.axhline(args.peptide_threshold, color='red', linestyle='--',
               label=f'Zscore Threshold-{args.peptide_threshold}')  # Add threshold line

    ax.set_xlabel("Samples", fontsize=15)  # Set x-axis label
    ax.set_ylabel("Z-Scores", fontsize=15)  # Set y-axis label
    plt.xticks(rotation=30, ha="right", fontsize=12)  # Rotate x-axis labels for better visibility
    plt.title("Boxplot of Enriched Peptides", fontsize=15)  # Set plot title
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1))  # Add legend to the plot

    output_file_path = os.path.join(args.output_path, f"{args.outputName}.pdf")  # Define output file path
    plt.savefig(output_file_path, dpi=300, bbox_inches="tight")  # Save the plot as a PDF

    if args.verbose:  # If verbose mode is enabled, print plot save status
        print(f"Plot saved to: {output_file_path}")

    plt.show()
