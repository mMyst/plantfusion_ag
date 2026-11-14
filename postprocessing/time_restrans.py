import numpy as np
import matplotlib.pyplot as plt
import os

# 1. Chargement de la tuile de base
filepath = r"C:\Users\agrumel\Code\Python_Ecophy\plantfusion_ag\outputs\WheatLuz_LUBBAC_emergence_adel\historique_restrans_caribu.npy"
filepath =  r"C:\Users\agrumel\Documents\Data\Simul LUBBAC\WheatLuz_LUBBAC_emergence_adel\historique_restrans_caribu.npy"

if not os.path.exists(filepath):
    raise FileNotFoundError(f"Fichier introuvable : {filepath}")

data = np.load(filepath)

surface_tuile_m2 = 0.09 
X_orig, Y_orig = data.shape[-2], data.shape[-1]
surface_case_m2 = surface_tuile_m2 / (X_orig * Y_orig)

data = data / surface_case_m2

# 2. Extraction de la plus grande couche Z (sommet de la canopée)
# Les dimensions sont (T, Z, X, Y). On prend le dernier index de Z.
data_top_z = data[:, -1, :, :] # Forme : (T_max, X_max, Y_max)

# 3. Subdivision en 7 colonnes
# On découpe l'axe Y (le dernier, axis=-1) en 7 parts.
# Index 0 à 6 correspondent aux colonnes 1 à 7.
colonnes = np.array_split(data_top_z, 7, axis=-1)

# 4. Regroupements spatiaux
# Luzerne : colonnes 1 et 7 (index 0 et 6)
luzerne = np.concatenate((colonnes[0], colonnes[6]), axis=-1)

# Blé : colonnes 4 et 5 (index 3 et 4)
ble = np.concatenate((colonnes[3], colonnes[4]), axis=-1)

# Colonnes restantes indépendantes
col2 = colonnes[1] # Colonne 2
col3 = colonnes[2] # Colonne 3
col6 = colonnes[5] # Colonne 6

# 5. Calcul des moyennes et écart-types (sur l'axe spatial X et sous-Y)
mean_luz = np.mean(luzerne, axis=(1, 2))
std_luz = np.std(luzerne, axis=(1, 2))

mean_ble = np.mean(ble, axis=(1, 2))
std_ble = np.std(ble, axis=(1, 2))

mean_c2 = np.mean(col2, axis=(1, 2))
std_c2 = np.std(col2, axis=(1, 2))

mean_c3 = np.mean(col3, axis=(1, 2))
std_c3 = np.std(col3, axis=(1, 2))

mean_c6 = np.mean(col6, axis=(1, 2))
std_c6 = np.std(col6, axis=(1, 2))

# 6. Création du graphique
temps = np.arange(data_top_z.shape[0])

plt.figure(figsize=(14, 8))

# Tracé Luzerne (Vert)
plt.plot(temps, mean_luz, label="Transmis Luzerne (Col 1 & 7)", color='#2ca02c', linewidth=2.5)
plt.fill_between(temps, mean_luz - std_luz, mean_luz + std_luz, color='#2ca02c', alpha=0.15)

# Tracé Blé (Orange/Or)
plt.plot(temps, mean_ble, label="Transmis Blé (Col 4 & 5)", color='#ff7f0e', linewidth=2.5)
plt.fill_between(temps, mean_ble - std_ble, mean_ble + std_ble, color='#ff7f0e', alpha=0.15)

# Tracés des colonnes restantes (tons de bleu/violet pour les distinguer)
plt.plot(temps, mean_c2, label="Transmis Col 2", color='#1f77b4', linestyle='--', linewidth=2)
plt.fill_between(temps, mean_c2 - std_c2, mean_c2 + std_c2, color='#1f77b4', alpha=0.1)

plt.plot(temps, mean_c3, label="Transmis Col 3", color='#9467bd', linestyle='--', linewidth=2)
plt.fill_between(temps, mean_c3 - std_c3, mean_c3 + std_c3, color='#9467bd', alpha=0.1)

plt.plot(temps, mean_c6, label="Transmis Col 6", color='#17becf', linestyle='--', linewidth=2)
plt.fill_between(temps, mean_c6 - std_c6, mean_c6 + std_c6, color='#17becf', alpha=0.1)

# Formatage
plt.title("Évolution de la transmission lumineuse par bande (Base 7 colonnes)", fontsize=15, pad=15)
plt.xlabel("Pas de temps (t)", fontsize=12)
plt.ylabel("Moyenne res_trans (W/m²)", fontsize=12) 
plt.legend(loc="upper right", fontsize=11, bbox_to_anchor=(1.25, 1))
plt.grid(True, linestyle=':', alpha=0.7)
plt.tight_layout()

# Affichage
plt.show()