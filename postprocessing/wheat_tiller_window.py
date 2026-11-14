import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# 0. Chargement du fichier de correspondance des dates
df_dates = pd.read_csv(r"C:\Users\agrumel\Code\Python_Ecophy\plantfusion_ag\inputs_fspmwheat\LUBBAC_H_24_25_capteurs.csv")
df_dates = pd.read_csv(r"C:/Users/agrumel/Code/Python_Ecophy/WheatFspm/fspm-wheat/example/Darwinkel/inputs/Darwinkel_open-meteo_since19101976.csv")


# Conversion de la colonne Date en objet datetime, puis extraction au format "jj/mm"
df_dates['Date_parsed'] = pd.to_datetime(df_dates['Date'], format='%d/%m/%Y', dayfirst=True)
df_dates['jj_mm'] = df_dates['Date_parsed'].dt.strftime('%d/%m')

# Création d'un dictionnaire rapide qui associe t -> "jj/mm"
t_to_date = dict(zip(df_dates['t'], df_dates['jj_mm']))

# 1. Dictionnaire des fichiers CSV à comparer : { "Nom simulation": "chemin/vers/fichier.csv" }
csv_files = {
    "Wheat": r"C:\Users\agrumel\Code\Python_Ecophy\plantfusion_ag\tiller_in_window_evolution_pure.csv",
    "WheatLuz": r"C:\Users\agrumel\Code\Python_Ecophy\plantfusion_ag\tiller_in_window_evolution_WL.csv",
    "latest": r"C:\Users\agrumel\Code\Python_Ecophy\plantfusion_ag\tiller_in_window_evolution.csv"
}

csv_files = {
   
    "Wheat": r"C:\Users\agrumel\Documents\Data\Simul LUBBAC\Wheat_LUBBAC_emergence_gai07\wheat\tiller_in_window_evolution.csv",
    "WheatLuzP" : r"C:\Users\agrumel\Documents\Data\Simul LUBBAC\WheatLuz_LUBBAC_emergence_typeK\wheat\tiller_in_window_evolution.csv",
    "WheatLuz": r"C:\Users\agrumel\Documents\Data\Simul LUBBAC\WheatLuz_LUBBAC_emergence_gai07\wheat\tiller_in_window_evolution.csv",
    "WheatTB" : r"C:\Users\agrumel\Documents\Data\Simul LUBBAC\WheatLuz_LUBBAC_emergence_TB\wheat\tiller_in_window_evolution.csv",
    
    #"Wheat": r"C:\Users\agrumel\Code\Python_Ecophy\plantfusion_ag\outputs\Wheat_LUBBAC_emergence_\wheat\tiller_in_window_evolution.csv",
    #"WheatLuz": r"C:\Users\agrumel\Code\Python_Ecophy\plantfusion_ag\outputs\WheatLuz_LUBBAC_emergence_gai07\wheat\tiller_in_window_evolution.csv",

}

csv_files = {
#"Dens. 100" : r"C:/Users/agrumel/Code/Python_Ecophy/WheatFspm/fspm-wheat/example/Darwinkel/darwinkel/GAIc_0.07_dens_100_delay_2.0_buf_1/outputs/tiller_in_window_evolution.csv",
"Dens. 200" : r"C:/Users/agrumel/Code/Python_Ecophy/WheatFspm/fspm-wheat/example/Darwinkel/darwinkel/GAIc_0.07_dens_200_delay_2.0_buf_1/outputs/tiller_in_window_evolution.csv",
"Dens. 250" : r"C:/Users/agrumel/Code/Python_Ecophy/WheatFspm/fspm-wheat/example/Darwinkel/darwinkel/GAIc_0.07_dens_250_delay_2.0_buf_1/outputs/tiller_in_window_evolution.csv",
"Dens. 400" : r"C:/Users/agrumel/Code/Python_Ecophy/WheatFspm/fspm-wheat/example/Darwinkel/darwinkel/GAIc_0.07_dens_400_delay_2.0_buf_1/outputs/tiller_in_window_evolution.csv",
"Dens. 800" : r"C:/Users/agrumel/Code/Python_Ecophy/WheatFspm/fspm-wheat/example/Darwinkel/darwinkel/GAIc_0.07_dens_800_delay_2.0_buf_1/outputs/tiller_in_window_evolution.csv",

}



# Chargement des données (séparateur virgule par défaut)
dfs = {name: pd.read_csv(path, sep=",") for name, path in csv_files.items()}

# Liste des talles
all_tillers = sorted(
    list(set(col for df in dfs.values() for col in df.columns if col.startswith("T"))),
    key=lambda x: int(x[1:]) if x[1:].isdigit() else x
)


# =======================================================
# 3. EXTRACTION DES DATES POUR LE TABLEAU
# =======================================================
summary_data = []

for sim_name, df in dfs.items():
    row_data = {"Simulation": sim_name}
    
    for tiller in all_tillers:
        if tiller in df.columns:
            # On récupère tous les 't' où la talle est in_window (True)
            active_t = df[df[tiller].fillna(False).astype(bool)]['t']
            
            if not active_t.empty:
                start_t = active_t.min()
                end_t = active_t.max()
                
                # Conversion du t en jj/mm via le dictionnaire
                start_date = t_to_date.get(int(start_t), f"t={int(start_t)}")
                end_date = t_to_date.get(int(end_t), f"t={int(end_t)}")
                
                row_data[tiller] = f"{start_date} au {end_date}"
            else:
                row_data[tiller] = "-"  # Jamais active
        else:
            row_data[tiller] = "-"
            
    summary_data.append(row_data)

# Création d'un DataFrame propre pour le tableau
df_summary = pd.DataFrame(summary_data)


# =======================================================
# 4. CRÉATION DU GRAPHIQUE (GRAPHIQUES + TABLEAU)
# =======================================================
colors = plt.get_cmap("tab10").colors

# Utilisation de GridSpec pour diviser la figure en 3 zones
fig = plt.figure(figsize=(11, 8.5)) # Figure un peu plus haute
gs = fig.add_gridspec(3, 1, height_ratios=[1, 2.5, 0.8]) # Le ratio 0.8 est pour le tableau

ax_gaip = fig.add_subplot(gs[0])
ax_tillers = fig.add_subplot(gs[1], sharex=ax_gaip)
ax_table = fig.add_subplot(gs[2]) # Axe réservé au tableau

num_sims = len(dfs)
offset = 0.8 / num_sims 

# -- TRACÉS DES COURBES ET DES BARRES --
for s_idx, (sim_name, df) in enumerate(dfs.items()):
    c = colors[s_idx % len(colors)]
    
    if "GAIp" in df.columns:
        ax_gaip.plot(df["t"], df["GAIp"], label=sim_name, color=c, lw=2)
    
    for t_idx, tiller in enumerate(all_tillers):
        if tiller not in df.columns: continue
            
        y_pos = t_idx + (s_idx - (num_sims - 1) / 2.0) * offset
        y_series = df[tiller].fillna(False).astype(bool).map({True: y_pos, False: np.nan})
        
        ax_tillers.plot(
            df["t"], y_series, color=c, linewidth=10, 
            solid_capstyle="butt", label=sim_name if t_idx == 0 else ""
        )

# -- FINITIONS DES AXES --
for t_idx in range(len(all_tillers) - 1):
    ax_tillers.axhline(t_idx + 0.5, color="gray", linestyle=":", alpha=0.3)

ax_tillers.set_yticks(range(len(all_tillers)))
ax_tillers.set_yticklabels(all_tillers)
ax_tillers.set_ylabel("Tillers")
ax_gaip.set_ylabel("GAIp")
ax_tillers.set_xlabel("Date (dd/mm)")

# Formateur de dates sur l'axe X
def format_date(x, pos):
    return t_to_date.get(int(x), "")

ax_tillers.xaxis.set_major_formatter(ticker.FuncFormatter(format_date))
plt.setp(ax_tillers.xaxis.get_majorticklabels(), rotation=45, ha='right')

ax_gaip.legend(loc="upper left")
ax_tillers.legend(loc="upper left")


# =======================================================
# 5. INTÉGRATION DU TABLEAU DANS LA FIGURE
# =======================================================
ax_table.axis('off') # On cache les axes (cadre, graduations) de ce subplot

# Création du tableau matplotlib
table = ax_table.table(
    cellText=df_summary.values,
    colLabels=df_summary.columns,
    cellLoc='center',
    loc='center',
    colColours=['#f2f2f2'] * len(df_summary.columns) # Fond gris clair pour l'en-tête
)

# Style du tableau
table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1, 1.6) # On étire légèrement le tableau pour qu'il soit plus aéré

plt.tight_layout()
plt.savefig("plot_with_table.png", dpi=300)
plt.show()