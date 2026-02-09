import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as mticker
from matplotlib.colors import PowerNorm

def display_root_profile(csv_file, layer_thickness=5):
    """
    Reads a CSV file, calculates the root length per soil layer over time,
    and displays the evolution as a heatmap.

    Args:
        csv_file (str): Path to the 'outputs_rootslog.csv' file.
        layer_thickness (float): Thickness of a soil layer in cm.
    """
    try:
        df = pd.read_csv(csv_file)
    except FileNotFoundError:
        print(f"Error: File not found at '{csv_file}'.")
        return
    except Exception as e:
        print(f"An error occurred while reading the CSV file: {e}")
        return

    # Check if required columns exist
    required_columns = ['t', 'x', 'y', 'z', 'roots_length']
    if not all(col in df.columns for col in required_columns):
        print(f"Error: CSV file must contain columns: {required_columns}")
        return

    # Group by time (t) and soil layer (z) and sum the root lengths
    root_profile = df.groupby(['t', 'z'])['roots_length'].sum().reset_index()

    root_profile = root_profile.pivot(index='t', columns='z', values='roots_length').fillna(0)

    # Calculate normalized profile
    root_profile_norm = root_profile.div(root_profile.sum(axis=1), axis=0).fillna(0)

    # Create heatmaps
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12), sharex=True)
    extent = [root_profile.index.min() - 0.5, root_profile.index.max() + 0.5,
              root_profile.columns.min() - 0.5, root_profile.columns.max() + 0.5]

    def layer_to_depth(x):
        return (x + 0.5) * layer_thickness

    def depth_to_layer(x):
        return (x / layer_thickness) - 0.5
    
    # Plot 1: Absolute Root Length
    im1 = ax1.imshow(root_profile.T, aspect='auto', origin='lower', cmap='viridis', extent=extent, norm=PowerNorm(gamma=0.3))  # Transpose the data to have z on y-axis
    ax1.set_ylabel('Soil Layer (z)')
    ax1.set_title('Root Length per Soil Layer over Time')
    fig.colorbar(im1, ax=ax1, label='Root Length')
    ax1.invert_yaxis()
    ax1.set_yticks(ticks=range(root_profile.columns.min(), root_profile.columns.max() + 1))
    
    secax1 = ax1.secondary_yaxis(-0.1, functions=(layer_to_depth, depth_to_layer))
    secax1.set_ylabel('Depth (cm)')

    # Plot 2: Normalized Root Length
    im2 = ax2.imshow(root_profile_norm.T, aspect='auto', origin='lower', cmap='viridis', extent=extent, norm=PowerNorm(gamma=0.3))
    ax2.set_xlabel('Time (DOY)')
    ax2.set_ylabel('Soil Layer (z)')
    ax2.set_title('Normalized Root Length per Soil Layer over Time')
    fig.colorbar(im2, ax=ax2, label='Normalized Root Length')
    ax2.invert_yaxis()
    ax2.set_yticks(ticks=range(root_profile.columns.min(), root_profile.columns.max() + 1))

    secax2 = ax2.secondary_yaxis(-0.1, functions=(layer_to_depth, depth_to_layer))
    secax2.set_ylabel('Depth (cm)')

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    # Example usage:


    #csv_file = r"C:\Users\agrumel\Code\Python_Ecophy\plantfusion_ag\outputs\cnwheat_soil3ds_debug\profile\0.2m\soil\outputs_rootslog.csv" 

    csv_file = r"C:\Users\agrumel\Code\Python_Ecophy\plantfusion_ag\outputs\cnwheat_soil3ds\homogeneous\1.5m\soil\outputs_rootslog.csv" 
    display_root_profile(csv_file, layer_thickness=5)
