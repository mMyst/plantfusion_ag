import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import warnings

# --- 1. Parameters, Configuration and Paths ---
CHEMIN_BASE = r"C:\Users\agrumel\Documents\Data\Simul LUBBAC"
CHEMIN_SUFFIXE = "historique_restrans_caribu.npy"
filepath_xls = r"inputs_soil_legume\LUBBAC_D_24_25_capteurs.xls"

# Paths to measured data
path_rayo_incid = r"//pnas1.stockage.inra.fr/nap-urp3f/Projets/2024_Lubbac_GLouarn/meteo/LUBBAC_Rayo_air_D.csv"
path_rayo_trans = r"//pnas1.stockage.inra.fr/nap-urp3f/Projets/2024_Lubbac_GLouarn/meteo/LUBBAC_Rayo_bacs_D.csv"
path_ref_incid  = r"//pnas1.stockage.inra.fr/nap-urp3f/Projets/2024_Lubbac_GLouarn/plan bacs capteurs.xlsx"
path_ref_lubbac = r"//pnas1.stockage.inra.fr/nap-urp3f/Projets/2024_Lubbac_GLouarn/ref_LUBBAC.xlsx"

# Mapping configs to columns and modalities
CONFIG_MAPPING = {
    "WheatLuz_LUBBAC_emergence_adel": {"col": 0, "modality": "LUZ", "title": "LUZ"},
    "WheatLuzREG_LUBBAC_emergence_adel": {"col": 1, "modality": "REG_LUZ", "title": "REG_LUZ"}
}

if not os.path.exists(filepath_xls):
    raise FileNotFoundError(f"File not found : {filepath_xls}")

# --- 2. Load Global Weather ---
df_meteo_global = pd.read_excel(filepath_xls)
df_meteo_global['Date'] = pd.to_datetime(
    df_meteo_global['year'].astype(str) + '-' + 
    df_meteo_global['month'].astype(str) + '-' + 
    df_meteo_global['day'].astype(str)
)

annee_simul = df_meteo_global['year'].iloc[0]
date_debut = pd.to_datetime(f"{annee_simul}-09-25")
df_meteo_global = df_meteo_global[df_meteo_global['Date'] >= date_debut].reset_index(drop=True)

# Conversion factor W/m² -> PPFD 
facteur_conversion = 4.57 

# Zeta Function
def calculer_zeta(tr_ppfd):
    a = 3.09
    b = 1.59
    x1 = 0
    x2 = 2
    c = 0
    d = 1.12
    tr_ppfd_clip = np.clip(tr_ppfd, x1, x2)
    numerateur = 1 - np.exp(-a * (tr_ppfd_clip - x1))
    denominateur = 1 - np.exp(-a * (x2 - x1))
    terme_interne = (c**b) + ((d**b) - (c**b)) * (numerateur / denominateur)
    return np.maximum(terme_interne, 0) ** (1/b)

# --- 3. Prepare Measured Data ---
# Load references
ref_incid = pd.read_excel(path_ref_incid, sheet_name='ref')
ref_incid['ID_sensor'] = ref_incid['ID_sensor'].astype(str)
ref_lubbac = pd.read_excel(path_ref_lubbac)
ref_lubbac['ID_BAC'] = ref_lubbac['ID_BAC'].astype(str)

# Process Transmitted Radiation (Bacs)
rayo_transm = pd.read_csv(path_rayo_trans, sep=';', low_memory=False)
noms_colonnes = rayo_transm.columns + "_" + rayo_transm.iloc[3].astype(str)
rayo_transm.columns = noms_colonnes
rayo_transm = rayo_transm.iloc[37:].reset_index(drop=True)

pivot_transm = rayo_transm.melt(id_vars=['Variable_Localisation 4'], var_name="sensor", value_name="rayo")
pivot_transm = pivot_transm.rename(columns={'Variable_Localisation 4': 'Date'})
pivot_transm['Date'] = pd.to_datetime(pivot_transm['Date']).dt.normalize()
pivot_transm['rayo'] = pd.to_numeric(pivot_transm['rayo'], errors='coerce')

# Extract IDs
pivot_transm[['ID_sensor', 'bac_id']] = pivot_transm['sensor'].str.split('_', n=1, expand=True)
pivot_transm['ID_sensor'] = pivot_transm['ID_sensor'].str.replace("L", "L-")
pivot_transm['bac_id'] = pivot_transm['bac_id'].str.replace("bac", "")

# Merge and Filter
pivot_transm = pivot_transm.merge(ref_incid, on="ID_sensor", how="left")
pivot_transm = pivot_transm.merge(ref_lubbac, left_on="bac_id", right_on="ID_BAC", how="left")

flo_date = pd.to_datetime("2025-05-15")
start_wheat = pd.to_datetime("2024-11-18")
pivot_transm = pivot_transm.dropna(subset=['rayo'])
pivot_transm = pivot_transm[(pivot_transm['Date'] != pd.to_datetime("2024-10-22")) & 
                            (pivot_transm['Date'] < flo_date) & 
                            (pivot_transm['Date'] > start_wheat)]

# Standardize Alfalfa naming
pivot_transm['Type 2'] = pivot_transm['Type 2'].replace({'Luz': 'Alf'})

# Process Incident Radiation (to get daily mean)
rayo_incid = pd.read_csv(path_rayo_incid, sep=';', low_memory=False)
rayo_incid = rayo_incid.iloc[74:].reset_index(drop=True)
rayo_incid = rayo_incid.melt(id_vars=['Variable'], var_name="sensor", value_name="Rayo")
rayo_incid = rayo_incid.rename(columns={'Variable': 'Date'})
rayo_incid['Date'] = pd.to_datetime(rayo_incid['Date']).dt.normalize()
rayo_incid['Rayo'] = pd.to_numeric(rayo_incid['Rayo'], errors='coerce')
rayo_incid.loc[rayo_incid['Rayo'] > 2000, 'Rayo'] = np.nan
rayo_incid = rayo_incid[rayo_incid['Date'] < flo_date]

rayo_inci_mean = rayo_incid.groupby('Date', as_index=False)['Rayo'].mean()
rayo_inci_mean = rayo_inci_mean.rename(columns={'Rayo': 'inci_rayo'})

# Merge individual sensors with incident rayo to calculate per-sensor ratios
pivot_transm = pivot_transm.merge(rayo_inci_mean, on='Date', how='inner')
pivot_transm['ratio_transm'] = pivot_transm['rayo'] / pivot_transm['inci_rayo']
pivot_transm['zeta_mesure'] = calculer_zeta(pivot_transm['ratio_transm'])

# Aggregate Mean and Std per Modality, Type 2, and Date
grouped_stats = pivot_transm.groupby(['Modality', 'Type 2', 'Date']).agg(
    rayo_mean=('rayo', 'mean'),
    rayo_std=('rayo', 'std'),
    ratio_mean=('ratio_transm', 'mean'),
    ratio_std=('ratio_transm', 'std'),
    zeta_mean=('zeta_mesure', 'mean'),
    zeta_std=('zeta_mesure', 'std')
).reset_index()

# Determine the very first measured date to adjust x-axis later
first_measured_date = grouped_stats['Date'].min()

# --- 4. Plotting Setup ---
# 3 rows (PPFD, Ratio, Zeta), 2 columns (LUZ, REG_LUZ)
fig, axs = plt.subplots(3, 2, figsize=(20, 15), sharex=True)
fig.subplots_adjust(right=0.85) # Leave space for external legends

colors = {'Alf': '#2ca02c', 'Wht': '#ff7f0e'}
labels = {'Alf': 'Alfalfa', 'Wht': 'Wheat'}

# Dictionaries to collect unique handles for row legends
handles_row0, labels_row0 = {}, {}
handles_row1, labels_row1 = {}, {}
handles_row2, labels_row2 = {}, {}

def update_legend(handles_dict, labels_dict, ax):
    h, l = ax.get_legend_handles_labels()
    for handle, label in zip(h, l):
        if label not in labels_dict.values():
            handles_dict[label] = handle
            labels_dict[label] = label

# --- 5. Loop over simulated configurations ---
for config_name, info in CONFIG_MAPPING.items():
    filepath_npy = os.path.join(CHEMIN_BASE, config_name, CHEMIN_SUFFIXE)
    col = info["col"]
    
    if not os.path.exists(filepath_npy):
        print(f"Warning: File not found for config {config_name} -> {filepath_npy}")
        continue
        
    # Load NPY
    data = np.load(filepath_npy)
    surface_tuile_m2 = 0.09 
    X_orig, Y_orig = data.shape[-2], data.shape[-1]
    surface_case_m2 = surface_tuile_m2 / (X_orig * Y_orig)
    data = data / surface_case_m2

    data_top_z = data[:, -1, :, :] 
    T_max = data_top_z.shape[0]

    # Slice weather
    df_meteo = df_meteo_global.iloc[:T_max]
    dates_plot = df_meteo['Date'].values
    i0_incident_w = df_meteo['I0'].values
    i0_incident_ppfd = i0_incident_w * facteur_conversion
    i0_safe = np.where(i0_incident_ppfd == 0, np.nan, i0_incident_ppfd)

    # Spatial processing
    colonnes = np.array_split(data_top_z, 7, axis=-1)
    
    # Concatenate specific columns in 3D (T, Y, X_subset)
    luzerne = np.concatenate((colonnes[0], colonnes[6]), axis=-1)
    ble = np.concatenate((colonnes[3], colonnes[4]), axis=-1)

    # Convert per-cell values to PPFD
    luzerne_ppfd = luzerne * facteur_conversion
    ble_ppfd = ble * facteur_conversion

    # 1. PPFD Means and STDs
    mean_luz_ppfd = np.mean(luzerne_ppfd, axis=(1, 2))
    std_luz_ppfd = np.std(luzerne_ppfd, axis=(1, 2))
    mean_ble_ppfd = np.mean(ble_ppfd, axis=(1, 2))
    std_ble_ppfd = np.std(ble_ppfd, axis=(1, 2))

    # Broadcast i0_safe to 3D array (T, 1, 1) to divide element-wise
    i0_safe_3d = i0_safe[:, np.newaxis, np.newaxis]

    # Calculate ratios per cell
    ratio_luz_array = luzerne_ppfd / i0_safe_3d
    ratio_ble_array = ble_ppfd / i0_safe_3d

    # Calculate zeta per cell
    zeta_luz_array = calculer_zeta(ratio_luz_array)
    zeta_ble_array = calculer_zeta(ratio_ble_array)

    # 2. Ratio & Zeta Means and STDs (ignoring NaNs when I0 = 0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        
        mean_ratio_luz = np.nanmean(ratio_luz_array, axis=(1, 2))
        std_ratio_luz = np.nanstd(ratio_luz_array, axis=(1, 2))
        mean_ratio_ble = np.nanmean(ratio_ble_array, axis=(1, 2))
        std_ratio_ble = np.nanstd(ratio_ble_array, axis=(1, 2))

        mean_zeta_luz = np.nanmean(zeta_luz_array, axis=(1, 2))
        std_zeta_luz = np.nanstd(zeta_luz_array, axis=(1, 2))
        mean_zeta_ble = np.nanmean(zeta_ble_array, axis=(1, 2))
        std_zeta_ble = np.nanstd(zeta_ble_array, axis=(1, 2))

    # Plot 1: PPFD
    axs[0, col].plot(dates_plot, i0_incident_ppfd, label="Incident I0 (PPFD)", color='black', linestyle=':', linewidth=1.5)
    
    axs[0, col].plot(dates_plot, mean_luz_ppfd, label="Simulated Alfalfa", color=colors['Alf'], linestyle='-', linewidth=2)
    axs[0, col].fill_between(dates_plot, np.maximum(0, mean_luz_ppfd - std_luz_ppfd), mean_luz_ppfd + std_luz_ppfd, color=colors['Alf'], alpha=0.2)
    
    axs[0, col].plot(dates_plot, mean_ble_ppfd, label="Simulated Wheat", color=colors['Wht'], linestyle='-', linewidth=2)
    axs[0, col].fill_between(dates_plot, np.maximum(0, mean_ble_ppfd - std_ble_ppfd), mean_ble_ppfd + std_ble_ppfd, color=colors['Wht'], alpha=0.2)
    update_legend(handles_row0, labels_row0, axs[0, col])

    # Plot 2: Ratio
    axs[1, col].plot(dates_plot, mean_ratio_luz, label="Simulated Alfalfa", color=colors['Alf'], linestyle='-', linewidth=2)
    axs[1, col].fill_between(dates_plot, np.maximum(0, mean_ratio_luz - std_ratio_luz), mean_ratio_luz + std_ratio_luz, color=colors['Alf'], alpha=0.2)
    
    axs[1, col].plot(dates_plot, mean_ratio_ble, label="Simulated Wheat", color=colors['Wht'], linestyle='-', linewidth=2)
    axs[1, col].fill_between(dates_plot, np.maximum(0, mean_ratio_ble - std_ratio_ble), mean_ratio_ble + std_ratio_ble, color=colors['Wht'], alpha=0.2)
    update_legend(handles_row1, labels_row1, axs[1, col])

    # Plot 3: Zeta
    axs[2, col].plot(dates_plot, mean_zeta_luz, label="Simulated Alfalfa", color=colors['Alf'], linestyle='-', linewidth=2)
    axs[2, col].fill_between(dates_plot, np.maximum(0, mean_zeta_luz - std_zeta_luz), mean_zeta_luz + std_zeta_luz, color=colors['Alf'], alpha=0.2)
    
    axs[2, col].plot(dates_plot, mean_zeta_ble, label="Simulated Wheat", color=colors['Wht'], linestyle='-', linewidth=2)
    axs[2, col].fill_between(dates_plot, np.maximum(0, mean_zeta_ble - std_zeta_ble), mean_zeta_ble + std_zeta_ble, color=colors['Wht'], alpha=0.2)
    update_legend(handles_row2, labels_row2, axs[2, col])


# --- 6. Add measured data points (Dashed lines + SD Ribbons) ---
for mod, info in CONFIG_MAPPING.items():
    col = info["col"]
    mod_name = info["modality"]
    
    for type_2 in ["Alf", "Wht"]:
        df_plot = grouped_stats[(grouped_stats['Modality'] == mod_name) & (grouped_stats['Type 2'] == type_2)].sort_values('Date')
        
        if df_plot.empty and type_2 == "Alf":
            df_plot = grouped_stats[(grouped_stats['Modality'] == mod_name) & (grouped_stats['Type 2'] == 'Luz')].sort_values('Date')
            
        if df_plot.empty:
            continue
            
        c = colors[type_2]
        label_mesure = f"Measured {labels[type_2]}"
        
        # Fill NaN standard deviations with 0 for continuous plotting
        rayo_std = df_plot['rayo_std'].fillna(0)
        ratio_std = df_plot['ratio_std'].fillna(0)
        zeta_std = df_plot['zeta_std'].fillna(0)
        
        # Plot 1 : Measured PPFD 
        axs[0, col].plot(df_plot['Date'], df_plot['rayo_mean'], color=c, linestyle='--', label=label_mesure, linewidth=2)
        axs[0, col].fill_between(df_plot['Date'], 
                                 np.maximum(0, df_plot['rayo_mean'] - rayo_std), 
                                 df_plot['rayo_mean'] + rayo_std, 
                                 color=c, alpha=0.2)
        update_legend(handles_row0, labels_row0, axs[0, col])
                       
        # Plot 2 : Measured Transmission Ratio
        axs[1, col].plot(df_plot['Date'], df_plot['ratio_mean'], color=c, linestyle='--', label=label_mesure, linewidth=2)
        axs[1, col].fill_between(df_plot['Date'], 
                                 np.maximum(0, df_plot['ratio_mean'] - ratio_std), 
                                 df_plot['ratio_mean'] + ratio_std, 
                                 color=c, alpha=0.2)
        update_legend(handles_row1, labels_row1, axs[1, col])
                       
        # Plot 3 : Measured Empirical Zeta
        axs[2, col].plot(df_plot['Date'], df_plot['zeta_mean'], color=c, linestyle='--', label=label_mesure, linewidth=2)
        axs[2, col].fill_between(df_plot['Date'], 
                                 np.maximum(0, df_plot['zeta_mean'] - zeta_std), 
                                 df_plot['zeta_mean'] + zeta_std, 
                                 color=c, alpha=0.2)
        update_legend(handles_row2, labels_row2, axs[2, col])

# --- 7. Final Formatting ---
for col in [0, 1]:
    modality_title = "LUZ Modality" if col == 0 else "REG_LUZ Modality"
    
    # Filter x-axis to start at the first measured date
    if not pd.isna(first_measured_date):
        axs[0, col].set_xlim(left=first_measured_date)
    
    # Row 0
    axs[0, col].set_title(f"Daily PPFD Evolution ({modality_title})", fontsize=13)
    axs[0, col].grid(True, linestyle=':', alpha=0.7)
    
    # Row 1
    axs[1, col].axhline(y=1, color='red', linestyle=':', alpha=0.5, label="Ratio = 1") 
    update_legend(handles_row1, labels_row1, axs[1, col])
    axs[1, col].set_title(f"Light Transmission Ratio ({modality_title})", fontsize=13)
    axs[1, col].set_ylim(0, 2.1) 
    axs[1, col].grid(True, linestyle=':', alpha=0.7)
    
    # Row 2
    axs[2, col].set_xlabel("Date", fontsize=12)
    axs[2, col].set_title(f"Zeta Parameter Evolution ({modality_title})", fontsize=13)
    axs[2, col].grid(True, linestyle=':', alpha=0.7)
    
    plt.setp(axs[2, col].xaxis.get_majorticklabels(), rotation=45, ha='right')

axs[0, 0].set_ylabel("PPFD", fontsize=11)
axs[1, 0].set_ylabel("Ratio (Transmitted / Incident)", fontsize=11)
axs[2, 0].set_ylabel("Zeta", fontsize=11)

# Apply external legends per row
fig.legend(handles_row0.values(), labels_row0.keys(), loc='center left', bbox_to_anchor=(0.86, 0.83), fontsize=10, title="Row 1 Legend")
fig.legend(handles_row1.values(), labels_row1.keys(), loc='center left', bbox_to_anchor=(0.86, 0.5), fontsize=10, title="Row 2 Legend")
fig.legend(handles_row2.values(), labels_row2.keys(), loc='center left', bbox_to_anchor=(0.86, 0.17), fontsize=10, title="Row 3 Legend")

plt.show()