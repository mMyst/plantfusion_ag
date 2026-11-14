import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import math

def plot_leaf1_custom_simulations(
    simulation_files: dict,
    target_column: str = "length",
    meteo_file: str = "inputs/meteo.csv",
    output_dir: str = "."
):
    """
    Analyse une liste explicite de fichiers 'elements_postprocessing.csv' 
    et trace l'évolution d'une variable spécifique (Feuille 1) par talle.
    Affichage en Facet par Simulation.
    """
    
    if not simulation_files:
        print("❌ Aucun fichier fourni dans le dictionnaire simulation_files.")
        return

    print(f"📂 Analyse de l'évolution de '{target_column}' pour {len(simulation_files)} simulations...")
    
    # --- 1. CHARGEMENT DU FICHIER MÉTÉO ---
    meteo_df = None
    if os.path.exists(meteo_file):
        try:
            meteo_df = pd.read_csv(meteo_file, comment='#')
            date_col = 'Date' if 'Date' in meteo_df.columns else ('date' if 'date' in meteo_df.columns else None)
            
            if not date_col:
                possible_cols = [col for col in meteo_df.columns if 'time' in col.lower() or 'date' in col.lower()]
                if possible_cols: date_col = possible_cols[0]
            
            if 't' in meteo_df.columns and date_col:
                meteo_df['date_formatted'] = pd.to_datetime(meteo_df[date_col]).dt.strftime('%d/%m')
                print(f"📅 Fichier météo chargé. Correspondance t/date active.")
            else:
                meteo_df = None
        except Exception as e:
            print(f"⚠️ Erreur météo ({e}). Mode fallback activé.")
            meteo_df = None


    # --- 2. LECTURE ET EXTRACTION DES DONNÉES ---
    all_dataframes = []
    
    for sim_name, file_path in simulation_files.items():
        if not os.path.exists(file_path):
            print(f"⚠️ Fichier introuvable ignoré : {file_path}")
            continue
            
        try:
            df = pd.read_csv(file_path)
            expected_cols = ['t', 'axis', 'metamer', 'organ', target_column]
            if not all(col in df.columns for col in expected_cols):
                print(f"⚠️ Colonne manquante dans {file_path}")
                continue
                
            df_tillers = df[(df['axis'] != 'MS') & (df['metamer'] == 1) & (df['organ'] == 'blade')].copy()
            df_tillers[target_column] = pd.to_numeric(df_tillers[target_column], errors='coerce')
            df_tillers = df_tillers.dropna(subset=[target_column])
            
            if df_tillers.empty: continue
                
            df_tillers['axis_clean'] = df_tillers['axis'].str.replace(r"b'|'", "", regex=True)
            df_grouped = df_tillers.groupby(['t', 'axis_clean'])[target_column].sum().reset_index()
            
            # Ajout du nom de la simulation pour le facet
            df_grouped['Simulation'] = sim_name
            all_dataframes.append(df_grouped)
            
        except Exception as e:
            print(f"❌ Erreur sur {file_path} : {e}")

    if not all_dataframes:
        print(f"❌ Aucune donnée valide trouvée.")
        return

    # --- 3. COMPILATION ---
    master_df = pd.concat(all_dataframes, ignore_index=True)
    
    unique_sims = list(simulation_files.keys())  # Conserve l'ordre du dictionnaire
    unique_axes = sorted(master_df['axis_clean'].unique())
    
    # --- 4. CONFIGURATION DE LA GRILLE MATPLOTLIB ---
    n_sims = len(unique_sims)
    ncols = min(n_sims, 3) 
    nrows = math.ceil(n_sims / ncols) if n_sims > 0 else 1

    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, 
                             figsize=(5 * ncols + 2, 4 * nrows + 1), 
                             sharex=True, sharey=True, squeeze=False, 
                             constrained_layout=True)
    
    axes_flat = axes.flatten()
    colors = plt.cm.tab10.colors 
    color_map = {axis: colors[i % len(colors)] for i, axis in enumerate(unique_axes)}

    # --- 5. REMPLISSAGE DES GRAPHIQUES ---
    for idx, sim_name in enumerate(unique_sims):
        ax = axes_flat[idx]
        subset_sim = master_df[master_df['Simulation'] == sim_name]
        
        for axis_name in unique_axes:
            grp = subset_sim[subset_sim['axis_clean'] == axis_name].sort_values('t')
            if not grp.empty:
                ax.plot(grp['t'], grp[target_column], 
                        color=color_map[axis_name], 
                        linestyle='-',
                        linewidth=2 if axis_name == 'T1' else 1.5, 
                        alpha=0.85)

        ax.set_title(f"{sim_name}", fontweight='bold', fontsize=12, pad=35)
        ax.grid(True, linestyle='--', alpha=0.6)
        
        if idx >= len(unique_sims) - ncols: 
            ax.set_xlabel("Temps (t, heures)", fontsize=10)
        if idx % ncols == 0: 
            ax.set_ylabel(f"{target_column} (Feuille 1)", fontsize=10)

    # Masquer les sous-graphes vides
    for idx in range(n_sims, len(axes_flat)):
        axes_flat[idx].set_visible(False)

    # --- 6. DEUXIÈME AXE X (DATES) ---
    for i, ax in enumerate(axes_flat):
        if not ax.get_visible(): continue
            
        ax2 = ax.twiny()
        ax2.set_xlim(ax.get_xlim())
        
        primary_ticks = ax.get_xticks()
        valid_ticks, tick_labels = [], []
        
        for t in primary_ticks:
            if ax.get_xlim()[0] <= t <= ax.get_xlim()[1]:
                if meteo_df is not None:
                    try:
                        closest_idx = (meteo_df['t'] - t).abs().idxmin()
                        tick_labels.append(meteo_df.loc[closest_idx, 'date_formatted'])
                        valid_ticks.append(t)
                    except: pass
                else:
                    try:
                        start_date = pd.to_datetime("1976-10-19")
                        date_val = start_date + pd.Timedelta(hours=float(t))
                        tick_labels.append(date_val.strftime('%d/%m'))
                        valid_ticks.append(t)
                    except: pass
                    
        ax2.set_xticks(valid_ticks)
        ax2.set_xticklabels(tick_labels, fontsize=9)
        if i < ncols: ax2.set_xlabel("Date (jj/mm)", fontsize=10, labelpad=10)

    # --- 7. LÉGENDE ---
    legend_elements = [Line2D([0], [0], color='w', label=r'$\bf{Axe \ (Talle)}$')]
    for axis_name, color in color_map.items():
        legend_elements.append(Line2D([0], [0], color=color, lw=2, label=f'{axis_name}'))

    fig.legend(handles=legend_elements, loc='outside right center', fontsize=11, frameon=True)
    fig.suptitle(f"Comparaison de '{target_column}' par Simulation", fontsize=15, fontweight='bold')

    # Sauvegarde
    os.makedirs(output_dir, exist_ok=True)
    safe_col_name = target_column.replace("/", "_").replace("\\", "_").upper()
    output_image = os.path.join(output_dir, f"COMPARAISON_{safe_col_name}.png")
    
    plt.savefig(output_image, dpi=300, bbox_inches='tight')
    print(f"💾 Graphique sauvegardé sous : {output_image}")
    plt.show()


if __name__ == "__main__":
    
    # === PARAMÈTRES À MODIFIER ===
    COLONNE_A_TESTER = "length" 
    FICHIER_METEO = "inputs_fspmwheat\LUBBAC_H_24_25_capteurs.csv"
    DOSSIER_SORTIE = "postprocessing\graphiques"
    
    # Définis ici tes simulations : {"Titre affiché": "chemin/vers/fichier.csv"}
    MES_SIMULATIONS = {
        "Wheat": r"C:\Users\agrumel\Documents\Data\Simul LUBBAC\Wheat_LUBBAC_emergence_adel\wheat\postprocessing\elements_postprocessing.csv",
        "WheatLuz" : r"C:\Users\agrumel\Documents\Data\Simul LUBBAC\WheatLuz_LUBBAC_emergence_adel\wheat\postprocessing\elements_postprocessing.csv",
        "WheatLuz Prost.": r"C:\Users\agrumel\Documents\Data\Simul LUBBAC\WheatLuz_LUBBAC_emergence_typeK\wheat\postprocessing\elements_postprocessing.csv",
        "WheatLuz AutoReg": r"C:\Users\agrumel\Documents\Data\Simul LUBBAC\WheatLuzREG_LUBBAC_emergence_adel_autoreg10\wheat\postprocessing\elements_postprocessing.csv",

    }
    # ==============================
    
    plot_leaf1_custom_simulations(
        simulation_files=MES_SIMULATIONS,
        target_column=COLONNE_A_TESTER,
        meteo_file=FICHIER_METEO,
        output_dir=DOSSIER_SORTIE
    )