import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import PowerNorm

def load_and_process_profiles(csv_file, file_label=None, divisor=1.0):
    """
    Reads a CSV file, divides root values by a coefficient, and computes 
    absolute and normalized root profiles per plant_name (or whole file).

    Args:
        csv_file (str): Path to the CSV file.
        file_label (str, optional): Label for the file.
        divisor (float): Coefficient to divide 'roots_length' by.

    Returns:
        list[dict]: A list of dicts containing {'title', 'prof', 'prof_norm', 'divisor'}.
    """
    if file_label is None:
        file_label = os.path.basename(os.path.dirname(csv_file)) or os.path.basename(csv_file)

    try:
        df = pd.read_csv(csv_file)
    except FileNotFoundError:
        print(f"Error: File not found at '{csv_file}'.")
        return []
    except Exception as e:
        print(f"Error reading file '{csv_file}': {e}")
        return []

    required_columns = ['t', 'x', 'y', 'z', 'roots_length']
    if not all(col in df.columns for col in required_columns):
        print(f"Error: '{csv_file}' must contain columns: {required_columns}")
        return []

    if divisor == 0:
        raise ValueError("Divisor coefficient cannot be zero.")

    # Apply scaling coefficient
    df['roots_length'] = df['roots_length'] / divisor

    results = []

    # Check whether 'plant_name' column exists
    if 'plant_name' in df.columns:
        plant_names = df['plant_name'].dropna().unique()
        for plant in plant_names:
            sub_df = df[df['plant_name'] == plant]
            
            prof = sub_df.groupby(['t', 'z'])['roots_length'].sum().reset_index()
            prof = prof.pivot(index='t', columns='z', values='roots_length').fillna(0)
            prof_norm = prof.div(prof.sum(axis=1), axis=0).fillna(0)
            
            title = f"{file_label}\nPlant: {plant}"
            if divisor != 1.0:
                title += f" (÷ {divisor})"

            results.append({
                'title': title,
                'prof': prof,
                'prof_norm': prof_norm
            })
    else:
        prof = df.groupby(['t', 'z'])['roots_length'].sum().reset_index()
        prof = prof.pivot(index='t', columns='z', values='roots_length').fillna(0)
        prof_norm = prof.div(prof.sum(axis=1), axis=0).fillna(0)
        
        title = file_label
        if divisor != 1.0:
            title += f" (÷ {divisor})"

        results.append({
            'title': title,
            'prof': prof,
            'prof_norm': prof_norm
        })

    return results


def compare_root_profiles(csv_files, labels=None, divisors=None, layer_thickness=5, share_cbar=True, gamma=0.3):
    """
    Displays root profiles side-by-side: 1 column per file and per plant_name.

    Args:
        csv_files (list[str]): List of CSV file paths.
        labels (list[str], optional): Custom labels for each CSV file.
        divisors (float | list[float] | dict, optional): Coefficient(s) to divide values by.
            - If float: applies the same divisor to all files.
            - If list: maps 1-to-1 to csv_files list.
            - If dict: maps file index or label/path to a divisor value.
        layer_thickness (float): Soil layer thickness in cm.
        share_cbar (bool): Synchronize color scales across all subplots.
        gamma (float): Gamma exponent for PowerNorm.
    """
    n_files = len(csv_files)
    if labels is None:
        labels = [None] * n_files

    # Normalize 'divisors' input to a list per file
    if divisors is None:
        divisor_list = [1.0] * n_files
    elif isinstance(divisors, (int, float)):
        divisor_list = [float(divisors)] * n_files
    elif isinstance(divisors, (list, tuple)):
        divisor_list = [float(d) for d in divisors]
        if len(divisor_list) != n_files:
            raise ValueError(f"Length of 'divisors' ({len(divisor_list)}) must match 'csv_files' ({n_files}).")
    elif isinstance(divisors, dict):
        divisor_list = []
        for i, path in enumerate(csv_files):
            # Check by index, then path, then custom label
            d = divisors.get(i, divisors.get(path, divisors.get(labels[i], 1.0)))
            divisor_list.append(float(d))
    else:
        raise TypeError("divisors must be None, a number, a list of numbers, or a dict.")

    all_series = []
    for path, label, div in zip(csv_files, labels, divisor_list):
        series_list = load_and_process_profiles(path, file_label=label, divisor=div)
        all_series.extend(series_list)

    n_cols = len(all_series)
    if n_cols == 0:
        print("No valid datasets found to display.")
        return

    # Compute global color scale limits
    if share_cbar:
        max_abs = max(item['prof'].values.max() for item in all_series)
        norm_abs = PowerNorm(gamma=gamma, vmin=0, vmax=max_abs) if max_abs > 0 else None
        norm_rel = PowerNorm(gamma=gamma, vmin=0, vmax=1.0)
    else:
        norm_abs = None
        norm_rel = None

    fig, axes = plt.subplots(
        nrows=2, 
        ncols=n_cols, 
        figsize=(4.5 * n_cols, 10), 
        sharey='row', 
        squeeze=False
    )

    def layer_to_depth(x):
        return (x + 0.5) * layer_thickness

    def depth_to_layer(x):
        return (x / layer_thickness) - 0.5

    for col_idx, item in enumerate(all_series):
        prof = item['prof']
        prof_norm = item['prof_norm']
        title = item['title']

        extent = [
            prof.index.min() - 0.5, prof.index.max() + 0.5,
            prof.columns.min() - 0.5, prof.columns.max() + 0.5
        ]

        cur_norm_abs = norm_abs if norm_abs else PowerNorm(gamma=gamma, vmin=0, vmax=prof.values.max())
        cur_norm_rel = norm_rel if norm_rel else PowerNorm(gamma=gamma, vmin=0, vmax=prof_norm.values.max())

        # --- Row 1: Absolute Root Length ---
        ax_abs = axes[0, col_idx]
        im1 = ax_abs.imshow(prof.T, aspect='auto', origin='lower', cmap='viridis', extent=extent, norm=cur_norm_abs)
        ax_abs.set_title(f"{title}\nAbsolute Root Length", fontsize=10)
        ax_abs.invert_yaxis()
        ax_abs.set_yticks(ticks=range(int(prof.columns.min()), int(prof.columns.max()) + 1))
        fig.colorbar(im1, ax=ax_abs, orientation='vertical', fraction=0.046, pad=0.04)

        if col_idx == 0:
            ax_abs.set_ylabel('Soil Layer (z)')
            secax1 = ax_abs.secondary_yaxis(-0.15, functions=(layer_to_depth, depth_to_layer))
            secax1.set_ylabel('Depth (cm)')

        # --- Row 2: Normalized Root Length ---
        ax_norm = axes[1, col_idx]
        im2 = ax_norm.imshow(prof_norm.T, aspect='auto', origin='lower', cmap='viridis', extent=extent, norm=cur_norm_rel)
        ax_norm.set_title("Normalized Root Length", fontsize=10)
        ax_norm.set_xlabel('Time (DOY)')
        ax_norm.invert_yaxis()
        ax_norm.set_yticks(ticks=range(int(prof.columns.min()), int(prof.columns.max()) + 1))
        fig.colorbar(im2, ax=ax_norm, orientation='vertical', fraction=0.046, pad=0.04)

        if col_idx == 0:
            ax_norm.set_ylabel('Soil Layer (z)')
            secax2 = ax_norm.secondary_yaxis(-0.15, functions=(layer_to_depth, depth_to_layer))
            secax2.set_ylabel('Depth (cm)')

    plt.tight_layout()
    plt.show()



if __name__ == "__main__":
    # List of CSV file paths to compare
    csv_files = [
        # r"C:/Users/agrumel/Documents/Data/Simul LUBBAC/Wheat_LUBBAC_emergence_adel/soil/outputs_rootslog.csv",
        # r"C:/Users/agrumel/Documents/Data/Simul LUBBAC/WheatLuz_LUBBAC_emergence_adel/soil/outputs_rootslog.csv",

        r"C:\Users\agrumel\Code\Python_Ecophy\plantfusion_ag\outputs\cnwheat_soil3ds_debug\homogeneous\1.5m\soil\outputs_rootslog.csv",
        r"C:\Users\agrumel\Code\Python_Ecophy\plantfusion_ag\outputs\cnwheat_soil3ds_debug\profile\0.2m\soil\outputs_rootslog.csv"
    ]

    # Custom labels for each column (optional)
    labels = [
        # "Pure Wheat",
        # "Wheat-Alf.",
        "Homogeneous",
        "Profile"
    ]

    # You can supply divisors as:
    # 1. A list per file: [1.0, 10.0]
    # 2. A dict: {"Emergence Adel": 1.0, "Homogeneous 1.5m": 2.5}
    # 3. A single float for all files: 100.0
    divisors = None #[2.0, 1.0]


    compare_root_profiles(
        csv_files=csv_files,
        labels=labels,
        divisors=divisors,
        layer_thickness=5,
        share_cbar=True  # Synchronizes color scale across simulations
    )