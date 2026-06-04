#this is a script to compare inputs meteo files from either inputs_fspmwheat 

#1 load n meteo files, n being the number of meteo files to compare. 
#2 for each column they have in common, plot them on the same graph.

#3 display the graphs in a pdf file

import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def compare_meteo_files(meteo_paths, output_pdf='meteo_comparison.pdf'):
    """
    Compares multiple meteo CSV files by plotting common columns.
    
    Args:
        meteo_paths (dict): Dictionary where keys are labels and values are file paths.
        output_pdf (str): Name of the output PDF file.
    """
    # 1. Load meteo files
    dataframes = {}
    for label, path in meteo_paths.items():
        if os.path.exists(path):
            if path.endswith('.csv'):
                df = pd.read_csv(path)
            elif path.endswith(('.xls', '.xlsx')):
                df = pd.read_excel(path)
            else:
                print(f"Unsupported file format for {path}")
                continue

            # Ensure 't' or index is consistent for plotting
            if 'DOY' in df.columns:
                df['DOY'] = df['DOY'] + df['DOY'].shift(1).where(df['DOY'].diff() < 0, 0).cumsum()
                df.set_index('DOY', inplace=True)
            dataframes[label] = df
        else:
            print(f"Warning: Path {path} does not exist.")

    if not dataframes:
        print("No data to compare.")
        return

    # 2. Identify common columns
    common_columns = set.intersection(*(set(df.columns) for df in dataframes.values()))
    # Remove non-numeric columns like 'Date' for plotting
    numeric_cols = [col for col in common_columns if col.lower() != 'date']

    # Calculate common DOY limits based on the shortest window
    doy_min = max(df.index.min() for df in dataframes.values())
    doy_max = min(df.index.max() for df in dataframes.values())

    # 3. Plot and save to PDF
    with PdfPages(output_pdf) as pdf:
        for col in sorted(numeric_cols):
                 
            # plot for regular values
            plt.figure(figsize=(10, 6))
            for label, df in dataframes.items():
                mask = (df.index >= doy_min) & (df.index <= doy_max)
                df_filtered = df.loc[mask]
                plt.plot(df_filtered.index, df_filtered[col], label=label, alpha=0.7)

            plt.title(f'Comparison of {col}')
            plt.xlabel('DOY')
            plt.ylabel(col)
            plt.legend()
            plt.grid(True, linestyle='--', alpha=0.6)
            plt.tight_layout()
            pdf.savefig()
            plt.close()


            # plot for cumulative sums
            plt.figure(figsize=(10, 6))
            for label, df in dataframes.items():
                mask = (df.index >= doy_min) & (df.index <= doy_max)
                df_filtered = df.loc[mask]
                plt.plot(df_filtered.index, df_filtered[col].cumsum(), label=f'{label} (cumulative)', alpha=0.8)

            
            plt.title(f'Comparison of {col}')
            plt.xlabel('DOY')
            plt.ylabel(col)
            plt.legend()
            plt.grid(True, linestyle='--', alpha=0.6)
            plt.tight_layout()
            
            pdf.savefig()
            plt.close()

    print(f"Comparison PDF saved as {output_pdf}")

if __name__ == "__main__":
    # Define the paths to the meteo files you want to compare
    cnwheat_meteo_files = {
        'Ljutovac2002': os.path.join('inputs_fspmwheat', 'meteo_Ljutovac2002.csv'),
        'LUBBAC': os.path.join('inputs_fspmwheat', 'LUBBAC_H_24_25.csv'),

        # Add other paths as needed
    }

    legume_meteo_files = {
        'legume_lusignan': os.path.join('inputs_soil_legume', 'meteo_exemple.xls'),
        'LUBBAC': os.path.join('inputs_soil_legume', 'LUBBAC_D_24_25.xls'),

        # Add other paths as needed
    }


    # Run the comparison
    compare_meteo_files(legume_meteo_files, output_pdf='legume_meteo_comparison.pdf')
    compare_meteo_files(cnwheat_meteo_files, output_pdf='cnwheat_meteo_comparison.pdf')

