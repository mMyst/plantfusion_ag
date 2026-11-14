import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import os

# 1. Configuration: List your .npy files and their display names here
# fichiers_a_comparer = {
#     "Caribu": r"C:\Users\agrumel\Code\Python_Ecophy\plantfusion_ag\outputs\legume_CvR\historique_restrans_caribu.npy",     # Example: "Display Name": "filename.npy"
#     "Default": r"C:\Users\agrumel\Code\Python_Ecophy\plantfusion_ag\outputs\legume_CvR\historique_restrans_default.npy",
#     "bugrid" : r"C:\Users\agrumel\Code\Python_Ecophy\plantfusion_ag\outputs\legume_caribu_buggrid\historique_restrans_caribu.npy"
#     # Add as many as you want: "Experiment 3": "file3.npy",
# }

fichiers_a_comparer = {
    # "PureWheat": r"C:\Users\agrumel\Code\Python_Ecophy\plantfusion_ag\outputs\Wheat_LUBBAC_emergence_wind_gai\historique_restrans_caribu.npy",     # Example: "Display Name": "filename.npy"
    "WheatLuz": r"C:\Users\agrumel\Documents\Data\Simul LUBBAC\WheatLuz_LUBBAC_emergence_adel\historique_restrans_caribu.npy",
    "WheatLuzREG": r"C:\Users\agrumel\Documents\Data\Simul LUBBAC\WheatLuzREG_LUBBAC_emergence_adel\historique_restrans_caribu.npy",
    "Wheat": r"C:\Users\agrumel\Documents\Data\Simul LUBBAC\Wheat_LUBBAC_emergence_adel\historique_restrans_caribu.npy",

    
    # "Legume" : r"C:\Users\agrumel\Code\Python_Ecophy\plantfusion_ag\outputs\legume_LUBBAC_testparam\historique_restrans_caribu.npy"
    # Add as many as you want: "Experiment 3": "file3.npy",
}



# 2. Load and validate the data
data_list = []
titles = []

print("Loading data...")
for name, filepath in fichiers_a_comparer.items():
    if not os.path.exists(filepath):
        print(f"⚠️ Warning: The file {filepath} was not found and will be skipped.")
        continue
    
    data = np.load(filepath)
    
    # --- NOUVEAU CODE : CONVERSION W -> W/m2 ---
    # La tuile fait 30cm x 30cm, soit 0.09 m² au total
    surface_tuile_m2 = 0.09 
    
    # On récupère les dimensions spatiales X et Y d'origine
    X_orig, Y_orig = data.shape[-2], data.shape[-1]
    
    # On calcule la surface d'une seule subdivision
    surface_case_m2 = surface_tuile_m2 / (X_orig * Y_orig)
    
    # Conversion : W divisé par des m² donne des W/m²
    data = data / surface_case_m2
    # -------------------------------------------
    
    # Tuilage 3x3 (9 cases) : on répète 1x pour le temps, 1x pour Z, 3x pour X, 3x pour Y
    data = np.tile(data, (1, 1, 3, 3))
    
    data_list.append(data)
    titles.append(name)

if not data_list:
    raise ValueError("No valid .npy files were loaded. Check your file paths.")

# Verify that all loaded arrays have the exact same shape
base_shape = data_list[0].shape
for d, name in zip(data_list, titles):
    if d.shape != base_shape:
        raise ValueError(f"Shape mismatch! {name} has shape {d.shape}, expected {base_shape}")

T_max, Z_max, X_max, Y_max = base_shape
N_files = len(data_list)
print(f"Successfully loaded {N_files} files. Dimensions: {T_max} time steps, {Z_max} layers. Spatial grid (tiled 3x3): {X_max}x{Y_max}")

# 3. Calculate GLOBAL min and max for a unified color scale
vmin = min([np.min(d) for d in data_list])
vmax = max([np.max(d) for d in data_list])

# 4. Prepare the Matplotlib Figure with dynamic subplots
# Width scales with the number of files (5 inches per file)
fig, axes = plt.subplots(1, N_files, figsize=(5 * N_files, 6))

# Ensure 'axes' is always a list, even if there is only 1 file
if N_files == 1:
    axes = [axes]

# Leave space at the bottom for the sliders
plt.subplots_adjust(bottom=0.25, right=0.9) 

t_init = 0
z_init = Z_max - 1

images = []

# Draw the initial heatmaps
for i, ax in enumerate(axes):
    # Plot the first frame
    im = ax.imshow(data_list[i][t_init, z_init, :, :], 
                   cmap='viridis', origin='lower', vmin=vmin, vmax=vmax)
    
    ax.set_title(f"{titles[i]}\nTemps : {t_init} | Couche : {z_init}")
    ax.set_xlabel("Axe X")
    if i == 0:
        ax.set_ylabel("Axe Y") # Only show Y label on the first plot to avoid clutter
        
    images.append(im)

# Add a single colorbar shared across all subplots
fig.colorbar(images[0], ax=axes, label="Valeur res_trans (W/m²)", shrink=0.8, location='right')


# 5. Create the slider axes [left, bottom, width, height]
ax_t = plt.axes([0.15, 0.1, 0.7, 0.03])
ax_z = plt.axes([0.15, 0.05, 0.7, 0.03])

# Create the sliders
slider_t = Slider(ax_t, 'Temps (t)', 0, T_max - 1, valinit=t_init, valstep=1)
slider_z = Slider(ax_z, 'Couche (z)', 0, Z_max - 1, valinit=z_init, valstep=1)

# 6. Update function
def update(val):
    t = int(slider_t.val)
    z = int(slider_z.val)
    
    # Update every heatmap in the grid
    for i, im in enumerate(images):
        im.set_data(data_list[i][t, z, :, :])
        axes[i].set_title(f"{titles[i]}\nTemps : {t} | Couche : {z}")
        
    fig.canvas.draw_idle()

# Connect the update function to the sliders
slider_t.on_changed(update)
slider_z.on_changed(update)

# Show the interactive UI
plt.show()
