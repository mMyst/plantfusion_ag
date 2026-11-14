import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

# --- 1. Paramètres et Configuration ---
CHEMIN_BASE = r"C:\Users\agrumel\Documents\Data\Simul LUBBAC"
CHEMIN_SUFFIXE = "historique_restrans_caribu.npy"
filepath_xls = "inputs_soil_legume\LUBBAC_D_24_25_capteurs.xls"

CONFIG = [
    "WheatLuz_LUBBAC_emergence_adel",
    "WheatLuzREG_LUBBAC_emergence_adel"
]

# Couleurs attribuées à chaque configuration (ajoutez-en si vous avez plus de 2 configs)
COULEURS = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']

if not os.path.exists(filepath_xls):
    raise FileNotFoundError(f"Fichier introuvable : {filepath_xls}")

# --- 2. Chargement Météo Globale ---
df_meteo_global = pd.read_excel(filepath_xls)
df_meteo_global['Date'] = pd.to_datetime(
    df_meteo_global['year'].astype(str) + '-' + 
    df_meteo_global['month'].astype(str) + '-' + 
    df_meteo_global['day'].astype(str)
)

annee_simul = df_meteo_global['year'].iloc[0]
date_debut = pd.to_datetime(f"{annee_simul}-09-25")
df_meteo_global = df_meteo_global[df_meteo_global['Date'] >= date_debut].reset_index(drop=True)

# Facteur de conversion W/m² -> PPFD 
facteur_conversion = 4.57 

# Fonction Zeta
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

# --- 3. Préparation des graphiques ---
fig, axs = plt.subplots(3, 1, figsize=(14, 15), sharex=True)
incident_plotte = False # Pour ne tracer l'incident (I0) qu'une seule fois

# --- 4. Boucle sur les configurations ---
for idx, config_name in enumerate(CONFIG):
    filepath_npy = os.path.join(CHEMIN_BASE, config_name, CHEMIN_SUFFIXE)
    
    if not os.path.exists(filepath_npy):
        print(f"Attention: Fichier introuvable pour la config {config_name} -> {filepath_npy}")
        continue
        
    print(f"Traitement de la configuration : {config_name}")
    
    # Chargement NPY
    data = np.load(filepath_npy)
    surface_tuile_m2 = 0.09 
    X_orig, Y_orig = data.shape[-2], data.shape[-1]
    surface_case_m2 = surface_tuile_m2 / (X_orig * Y_orig)
    data = data / surface_case_m2

    data_top_z = data[:, -1, :, :] 
    T_max = data_top_z.shape[0]

    # Découpage de la météo pour la taille de cette simulation spécifique
    df_meteo = df_meteo_global.iloc[:T_max]
    dates_plot = df_meteo['Date'].values
    i0_incident_w = df_meteo['I0'].values
    i0_incident_ppfd = i0_incident_w * facteur_conversion
    i0_safe = np.where(i0_incident_ppfd == 0, np.nan, i0_incident_ppfd)

    # Traitement spatial
    colonnes = np.array_split(data_top_z, 7, axis=-1)
    luzerne = np.concatenate((colonnes[0], colonnes[6]), axis=-1)
    ble = np.concatenate((colonnes[3], colonnes[3]), axis=-1)

    mean_luz_w = np.mean(luzerne, axis=(1, 2))
    mean_ble_w = np.mean(ble, axis=(1, 2))

    mean_luz_ppfd = mean_luz_w * facteur_conversion
    mean_ble_ppfd = mean_ble_w * facteur_conversion

    ratio_luz = mean_luz_ppfd / i0_safe
    ratio_ble = mean_ble_ppfd / i0_safe

    zeta_luz = calculer_zeta(ratio_luz)
    zeta_ble = calculer_zeta(ratio_ble)
    
    # Style de tracé
    couleur = COULEURS[idx % len(COULEURS)]
    style_luzerne = '-'   # Trait plein
    style_ble = '--'      # Pointillés
    
    # Tracé 1 : PPFD
    if not incident_plotte:
        axs[0].plot(dates_plot, i0_incident_ppfd, label="I0 Incident (PPFD)", color='black', linestyle=':', linewidth=1.5)
        incident_plotte = True
        
    axs[0].plot(dates_plot, mean_luz_ppfd, label=f"Luzerne ({config_name})", color=couleur, linestyle=style_luzerne, linewidth=2)
    axs[0].plot(dates_plot, mean_ble_ppfd, label=f"Blé ({config_name})", color=couleur, linestyle=style_ble, linewidth=2)

    # Tracé 2 : Ratio
    axs[1].plot(dates_plot, ratio_luz, label=f"Luzerne ({config_name})", color=couleur, linestyle=style_luzerne, linewidth=2)
    axs[1].plot(dates_plot, ratio_ble, label=f"Blé ({config_name})", color=couleur, linestyle=style_ble, linewidth=2)

    # Tracé 3 : Zeta
    axs[2].plot(dates_plot, zeta_luz, label=f"Luzerne ({config_name})", color=couleur, linestyle=style_luzerne, linewidth=2)
    axs[2].plot(dates_plot, zeta_ble, label=f"Blé ({config_name})", color=couleur, linestyle=style_ble, linewidth=2)


# --- 5. Mise en forme finale des graphiques ---
# Ligne repère Ratio = 1
axs[1].axhline(y=1, color='red', linestyle=':', alpha=0.5, label="Ratio = 1") 

# Graphique 1
axs[0].set_ylabel("PPFD", fontsize=11)
axs[0].set_title("Évolution du PPFD Journalier (Transmis vs Incident)", fontsize=13)
axs[0].legend(loc="upper right", fontsize=9, bbox_to_anchor=(1.0, 1.0))
axs[0].grid(True, linestyle=':', alpha=0.7)

# Graphique 2
axs[1].set_ylabel("Ratio (Transmis / Incident)", fontsize=11)
axs[1].set_title("Ratio de transmission lumineuse (TrPPFD)", fontsize=13)
axs[1].set_ylim(0, 2.1) 
axs[1].legend(loc="upper right", fontsize=9)
axs[1].grid(True, linestyle=':', alpha=0.7)

# Graphique 3
axs[2].set_xlabel("Date", fontsize=12)
axs[2].set_ylabel("Zeta", fontsize=11)
axs[2].set_title("Évolution du paramètre Zeta en fonction du ratio de transmission", fontsize=13)
axs[2].legend(loc="upper right", fontsize=9)
axs[2].grid(True, linestyle=':', alpha=0.7)

# Rotation des dates
plt.setp(axs[2].xaxis.get_majorticklabels(), rotation=45, ha='right')

plt.tight_layout()
plt.show()