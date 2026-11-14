import os
import pandas as pd
import matplotlib.pyplot as plt


def plot_leaf1_plant_comparative(files_dict, 
                                 target_column="length",
                                 epsilon=0.0,
                                 meteo_file="inputs_fspmwheat/meteo_Ljutovac2002_since_sowing_on_19981015.csv",
                                 sowing_date=None,
                                 output_image="COMPARATIVE_PLANT_COUNT.png"):
    """
    Analyse un dictionnaire de fichiers 'elements_postprocessing.csv' pour tracer l'évolution 
    du NOMBRE TOTAL de feuilles 1 (metamer=1, organ='blade') actives à l'échelle de la plante.
    L'axe des abscisses affiche des incréments tous les 10 DAS au format jj/mm (DAS).
    Corrige les valeurs de 't' en fonction de la date de semis fournie ou déduite.
    """
    
    if not files_dict:
        print("❌ Le dictionnaire de fichiers est vide.")
        return

    print(f"📂 Analyse du NOMBRE TOTAL de feuilles 1 ({target_column} > {epsilon}) par plante...")
    
    # --- 1. CHARGEMENT DU FICHIER MÉTÉO ---
    meteo_df = None
    if os.path.exists(meteo_file):
        try:
            meteo_df = pd.read_csv(meteo_file, comment='#')
            
            date_col = None
            if 'Date' in meteo_df.columns:
                date_col = 'Date'
            elif 'date' in meteo_df.columns:
                date_col = 'date'
            else:
                possible_cols = [col for col in meteo_df.columns if 'time' in col.lower() or 'date' in col.lower()]
                if possible_cols:
                    date_col = possible_cols[0]
            
            if 't' in meteo_df.columns and date_col:
                date_series = meteo_df[date_col].astype(str)
                meteo_df['date_parsed'] = pd.to_datetime(date_series, format='%d/%m/%Y', errors='coerce')
                if meteo_df['date_parsed'].isna().all():
                    meteo_df['date_parsed'] = pd.to_datetime(date_series, errors='coerce')

                meteo_df['date_formatted'] = meteo_df['date_parsed'].dt.strftime('%d/%m')
                
                if not meteo_df['date_parsed'].isna().all():
                    print(f"📅 Fichier météo chargé. Correspondance t/date active via '{date_col}'.")
                else:
                    meteo_df = None
        except Exception as e:
            print(f"⚠️ Erreur de lecture météo ({e}).")
            meteo_df = None

    # --- 1bis. DÉTERMINATION DE LA DATE DE SEMIS ET t_sowing ---
    t_sowing = 0
    if sowing_date is not None:
        sowing_date_obj = pd.to_datetime(sowing_date)
        print(f"🌱 Date de semis (fournie manuellement) : {sowing_date_obj.strftime('%d/%m')}")
        if meteo_df is not None:
            match = meteo_df[meteo_df['date_parsed'] == sowing_date_obj]
            if not match.empty:
                t_sowing = match['t'].iloc[0]
            else:
                print("⚠️ Date de semis introuvable dans la météo. 't' non corrigé (t_sowing=0).")
    elif meteo_df is not None:
        sowing_date_obj = meteo_df['date_parsed'].min()
        match = meteo_df[meteo_df['date_parsed'] == sowing_date_obj]
        if not match.empty:
            t_sowing = match['t'].iloc[0]
        print(f"🌱 Date de semis (déduite de la météo) : {sowing_date_obj.strftime('%d/%m')}")
    else:
        sowing_date_obj = pd.to_datetime("1998-10-15")
        print(f"⚠️ Fallback de sécurité : {sowing_date_obj.strftime('%d/%m')}")

    print(f"🔧 Correction du temps appliquée : t = t - {t_sowing} heures.")

    # --- 2. LECTURE ET EXTRACTION POUR CHAQUE FICHIER ---
    all_dataframes = []

    for label, file_path in files_dict.items():
        if not os.path.exists(file_path):
            print(f"⚠️ Fichier introuvable : {file_path}")
            continue
            
        try:
            df = pd.read_csv(file_path)
            expected_cols = ['t', 'axis', 'metamer', 'organ', target_column]
            if not all(col in df.columns for col in expected_cols):
                continue
                
            df_filtered = df[(df['axis'] != 'MS') & (df['metamer'] == 1) & (df['organ'] == 'blade')].copy()
            df_filtered[target_column] = pd.to_numeric(df_filtered[target_column], errors='coerce')
            df_filtered = df_filtered.dropna(subset=[target_column])
            
            if df_filtered.empty:
                continue
                
            df_filtered['is_active'] = (df_filtered[target_column] > epsilon).astype(int)
            df_grouped = df_filtered.groupby(['t'])['is_active'].sum().reset_index()
            
            # --- CORRECTION DU 't' ---
            df_grouped['t'] = df_grouped['t'] - t_sowing
            # Optionnel: on peut filtrer pour ne garder que la période post-semis (t >= 0)
            df_grouped = df_grouped[df_grouped['t'] >= 0]
            
            df_grouped['Label'] = label  
            all_dataframes.append(df_grouped)
            
        except Exception as e:
            print(f"❌ Erreur lors de la lecture de {file_path} : {e}")

    if not all_dataframes:
        print(f"❌ Aucune données valides trouvées après correction.")
        return

    # --- 3. CRÉATION DU GRAPHIQUE ---
    master_df = pd.concat(all_dataframes, ignore_index=True)
    fig, ax = plt.subplots(figsize=(10, 6), constrained_layout=True)

    colors = plt.cm.tab10.colors 
    color_map = {label: colors[i % len(colors)] for i, label in enumerate(files_dict.keys())}

    for label in files_dict.keys():
        grp = master_df[master_df['Label'] == label].sort_values('t')
        if not grp.empty:
            ax.plot(grp['t'], grp['is_active'], color=color_map[label], linewidth=2.5, alpha=0.85, label=label)

    ax.set_title(f"Nombre total de talles actives par plante ({target_column} > {epsilon})", fontsize=14, fontweight='bold', pad=15)
    ax.set_ylabel("Limbes actifs par plante", fontsize=11)
    ax.grid(True, linestyle='--', alpha=0.6)

    # --- 4. AXE X : INCRÉMENTS TOUS LES 10 DAS ---
    t_min = master_df['t'].min()
    t_max = master_df['t'].max()
    
    das_min = int(t_min // 24)
    das_max = int(t_max // 24)
    
    start_das = (das_min // 10) * 10
    end_das = ((das_max // 10) + 1) * 10
    
    das_ticks = range(start_das, end_das + 1, 10)
    
    valid_ticks = []
    tick_labels = []
    
    for das in das_ticks:
        t_val = das * 24  # Correspond à la valeur corrigée (heures depuis le semis)
        
        if t_min <= t_val <= t_max:
            try:
                # Puisque t_val est déjà corrigé (t=0 au semis), on ajoute directement à sowing_date_obj
                date_val = sowing_date_obj + pd.Timedelta(hours=t_val)
                date_str = date_val.strftime('%d/%m')
            except:
                date_str = ""
                    
            valid_ticks.append(t_val)
            tick_labels.append(f"{date_str}\n(DAS {das})")
                
    ax.set_xticks(valid_ticks)
    ax.set_xticklabels(tick_labels, fontsize=10, rotation=0, ha='center')
    ax.set_xlim(max(0, t_min), t_max) # On assure que le graphe commence au semis (t=0) ou après
    ax.set_xlabel("Date (jj/mm) & DAS", fontsize=11, labelpad=10)

    # --- 5. LÉGENDE GLOBALE ---
    ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1), fontsize=11, frameon=True)

    plt.savefig(output_image, dpi=300, bbox_inches='tight')
    print(f"💾 Graphique sauvegardé sous : {output_image}")
    plt.show()


if __name__ == "__main__":
    
    COLONNE_A_TESTER = "length"  
    EPSILON = 0.01  
    
    fichiers_a_comparer =  {
#   "CN_Wheat_Asso" : "C:/Users/agrumel/Documents/Données/Simul LUBBAC/WheatLuz_LUBBAC_emergence/wheat/postprocessing/elements_postprocessing.csv"
  #,"CN_Wheat_AssoReg" = "C:/Users/agrumel/Documents/Données/Simul LUBBAC/WheatLuzREG_LUBBAC_emergence/wheat/postprocessing/elements_postprocessing.csv"
  "CN_Wheat_Pure" : "C:/Users/agrumel/Documents/Données/Simul LUBBAC/Wheat_LUBBAC_emergence/wheat/postprocessing/elements_postprocessing.csv"
  ,"Wind"  :"C:/Users/agrumel/Documents/Données/Simul LUBBAC/Wheat_LUBBAC_emergence_wind/wheat/postprocessing/elements_postprocessing.csv"
#   ,"ZWind" : "C:/Users/agrumel/Documents/Données/Simul LUBBAC/Wheat_LUBBAC_emergence_Zwind/wheat/postprocessing/elements_postprocessing.csv"
  ,"WheatLuz LUBBAC test": "outputs\WheatLuz_LUBBAC_emergence_wind_test_gai\wheat\postprocessing\elements_postprocessing.csv"

  }

    # fichiers_a_comparer = {
    #     "Darwinkel" : "C:/Users/agrumel/Code/Python_Ecophy/WheatFspm/fspm-wheat/example/Darwinkel/darwinkel/GAIc_0.07_dens_250_delay_2.0_buf_1/postprocessing/elements_postprocessing.csv"
    # }

    # fichiers_a_comparer = {
    #     "cnwheat_default" : "c:/Users/agrumel/Documents/Données/Sorties CNWheat/cnwheat_default_emergence/wheat/postprocessing/elements_postprocessing.csv"
    # }

    
    
    plot_leaf1_plant_comparative(
        files_dict=fichiers_a_comparer,
        target_column=COLONNE_A_TESTER,
        epsilon=EPSILON,
        meteo_file="inputs_fspmwheat\LUBBAC_H_24_25_capteurs.csv",
        sowing_date="2024-11-18",  # Si None, le script prendra la 1ère date du fichier météo
        output_image="Comparatif_talles_F1_lubbac.png"
    )

        
