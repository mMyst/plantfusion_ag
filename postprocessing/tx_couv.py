import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ==========================================
# 1. PARAMÈTRES ET CONFIGURATION
# ==========================================
CHEMIN_BASE = "C:/Users/agrumel/Documents/Data/Simul LUBBAC"
CHEMIN_SUFFIXE = "light/vtk"

CONFIG_ROUTAGE = {
    "WheatLuz_LUBBAC_emergence_adel": [2, 3, 4, 5, 6, 7, 8, 9],
    "Wheat_LUBBAC_emergence_adel": [1],
}

# Mise à jour avec le nouveau fichier
FICHIER_CSV = "resultats_taux_couv_detaille.csv"

DATE_SEMIS = pd.to_datetime("2024-11-18") # Jour 0 (DAS 0)
DATE_LIMITE_RMSE = "2025-03-15"  

X_DOMAINE = 0.3      
RESOLUTION = 0.001   
HAUTEUR_CIBLE = 600  
LARGEUR_CIBLE = 1150 
DECALAGE_X = 50      

DOSSIER_SAUVEGARDE = r"C:\Users\agrumel\Documents\Data\Simul LUBBAC\Sauvergardes_Figures\latest\Taux_couverture_par_espece"


# ==========================================
# 2. CALCUL DU TAUX DE COUVERTURE DEPUIS LE VTK
# ==========================================
def generer_vtk_taux_couverture(dossier_vtk, doy, X_dom, res, h_cible, w_cible, simulation_name):
    """
    Lit le fichier VTK et retourne le taux de couverture global et par espèce.
    """
    if not os.path.exists(dossier_vtk):
        return None

    fichiers = os.listdir(dossier_vtk)
    chemin_vtk = None
    doy_str = str(doy).zfill(3) 
    
    motif = rf"^.*_\d+{doy_str}13\.vtk$" 
    for f in fichiers:
        if re.match(motif, f):
            chemin_vtk = os.path.join(dossier_vtk, f)
            break
            
    if not chemin_vtk:
        return None
        
    with open(chemin_vtk, 'r') as f:
        texte = f.read()
    
    try:
        deb_pts = texte.find("POINTS")
        fin_pts = texte.find("CELLS")
        pts = np.array([float(x) for l in texte[deb_pts:fin_pts].strip().split('\n')[1:] for x in l.split()]).reshape(-1, 3)
        
        deb_cells = texte.find("CELLS")
        fin_cells = texte.find("CELL_TYPES")
        cells = np.array([int(x) for l in texte[deb_cells:fin_cells].strip().split('\n')[1:] for x in l.split()]).reshape(-1, 4)[:, 1:]
        
        # Extraction de l'espèce (VegetationType)
        veg_types = np.zeros(len(cells), dtype=int)
        deb_veg = texte.find("VegetationType")
        if deb_veg != -1:
            veg_data_str = texte[deb_veg:].split('\n', 1)[1]
            words = veg_data_str.split()
            veg_types = np.array([int(float(w)) for w in words[:len(cells)]])
            
    except Exception:
        return None

    taille_base = int(X_dom / res)
    z_buffer_base = np.full((taille_base, taille_base), -np.inf)
    species_buffer_base = np.full((taille_base, taille_base), -1, dtype=int) 
    
    def f_bary(v0, v1, p):
        return (v1[0] - v0[0]) * (p[1] - v0[1]) - (v1[1] - v0[1]) * (p[0] - v0[0])

    for dx in [-X_dom, 0.0, X_dom]:
        for dy in [-X_dom, 0.0, X_dom]:
            pts_dec = pts.copy()
            pts_dec[:, 0] += dx
            pts_dec[:, 1] += dy
            
            for tri_idx, tri in enumerate(cells):
                v0, v1, v2 = pts_dec[tri[0]], pts_dec[tri[1]], pts_dec[tri[2]]
                v0_p, v1_p, v2_p = v0[:2]/res, v1[:2]/res, v2[:2]/res
                v_type = veg_types[tri_idx]
                
                min_x = max(0, int(np.floor(min(v0_p[0], v1_p[0], v2_p[0]))))
                max_x = min(taille_base - 1, int(np.ceil(max(v0_p[0], v1_p[0], v2_p[0]))))
                min_y = max(0, int(np.floor(min(v0_p[1], v1_p[1], v2_p[1]))))
                max_y = min(taille_base - 1, int(np.ceil(max(v0_p[1], v1_p[1], v2_p[1]))))
                
                if min_x > max_x or min_y > max_y: continue
                aire = f_bary(v0_p, v1_p, v2_p)
                if aire == 0: continue
                    
                for py in range(min_y, max_y + 1):
                    for px in range(min_x, max_x + 1):
                        p = [px + 0.5, py + 0.5]
                        w0, w1, w2 = f_bary(v1_p, v2_p, p), f_bary(v2_p, v0_p, p), f_bary(v0_p, v1_p, p)
                        
                        est_dedans = (w0 >= 0 and w1 >= 0 and w2 >= 0) if aire > 0 else (w0 <= 0 and w1 <= 0 and w2 <= 0)
                        if est_dedans:
                            z = (w0 * v0[2] + w1 * v1[2] + w2 * v2[2]) / aire
                            if z > z_buffer_base[py, px]:
                                z_buffer_base[py, px] = z
                                species_buffer_base[py, px] = v_type

    species_buffer_base = np.roll(species_buffer_base, shift=DECALAGE_X, axis=1)
    
    reps_y = int(np.ceil(h_cible / taille_base))
    reps_x = int(np.ceil(w_cible / taille_base))
    
    img_tuilee_species = np.tile(species_buffer_base, (reps_y, reps_x))
    img_tuilee_species = img_tuilee_species[:h_cible, :w_cible] 
    
    # Cartographie des ID d'espèces
    is_wheatluz = "WheatLuz" in simulation_name
    mask_wheat = (img_tuilee_species == (1 if is_wheatluz else 0))
    mask_alfalfa = (img_tuilee_species == 0) if is_wheatluz else np.zeros_like(img_tuilee_species, dtype=bool)
    
    return {
        "Total": np.mean(img_tuilee_species != -1),
        "Wheat": np.mean(mask_wheat),
        "Alfalfa": np.mean(mask_alfalfa)
    }


# ==========================================
# 3. EXTRACTION DU CSV ET COMPARAISON
# ==========================================
if __name__ == "__main__":
    print("=======================================")
    print("COMPARAISON DÉTAILLÉE PAR ESPÈCE (CSV vs VTK)")
    print("=======================================")
    
    if not os.path.exists(DOSSIER_SAUVEGARDE):
        os.makedirs(DOSSIER_SAUVEGARDE)

    try:
        df_csv = pd.read_csv(FICHIER_CSV, sep=',')
        df_csv['Date_dt'] = pd.to_datetime(df_csv['Date'])
        df_csv['DAS'] = (df_csv['Date_dt'] - DATE_SEMIS).dt.days
        df_csv['DOY'] = df_csv['Date_dt'].dt.dayofyear
        df_csv['Date_Affiche'] = df_csv['Date_dt'].dt.strftime('%d/%m')
        
        print(f"✅ CSV détaillé chargé avec succès. {df_csv['DAS'].nunique()} jours d'observation détectés.")
    except Exception as e:
        print(f"🛑 Erreur fatale : Impossible de traiter '{FICHIER_CSV}' ({e})")
        exit()

    resultats = []

    for simulation, bacs_reels in CONFIG_ROUTAGE.items():
        dossier_vtk = os.path.join(CHEMIN_BASE, simulation, CHEMIN_SUFFIXE).replace('\\', '/')
        is_wheatluz = "WheatLuz" in simulation
        
        if not os.path.exists(dossier_vtk):
            continue
            
        print(f"\n>> Traitement : {simulation} (Bacs: {bacs_reels})")
        df_simul = df_csv[df_csv['Bac'].isin(bacs_reels)]
        
        if df_simul.empty: continue
            
        # Agrégation des nouvelles colonnes
        df_obs_stats = df_simul.groupby(['DAS', 'DOY', 'Date_dt', 'Date_Affiche']).agg(
            Obs_Total_moy=('Moyenne_Total', 'mean'), Obs_Total_std=('SD_Total', 'mean'),
            Obs_Side_moy=('Moyenne_Side_Carre', 'mean'), Obs_Side_std=('SD_Side_Carre', 'mean'),
            Obs_Center_moy=('Moyenne_Center_Carre', 'mean'), Obs_Center_std=('SD_Center_Carre', 'mean'),
            Obs_IR_moy=('Moyenne_IR_Carre', 'mean'), Obs_IR_std=('SD_IR_Carre', 'mean')
        ).reset_index().sort_values(by='DAS')

        for _, row in df_obs_stats.iterrows():
            taux_simul = generer_vtk_taux_couverture(dossier_vtk, int(row['DOY']), X_DOMAINE, RESOLUTION, HAUTEUR_CIBLE, LARGEUR_CIBLE, simulation)
            
            if taux_simul is not None:
                diff_tot = taux_simul["Total"] - row['Obs_Total_moy']
                print(f"  - DAS {int(row['DAS']):3d} | Sim Total: {taux_simul['Total']:.3f} | Obs Total: {row['Obs_Total_moy']:.3f} | Diff: {diff_tot:+.3f}")
                
                resultats.append({
                    "Simulation": simulation,
                    "DAS": int(row['DAS']),
                    "Date_dt": row['Date_dt'],
                    "Date": row['Date_Affiche'],
                    
                    "Sim_Total": taux_simul["Total"],
                    "Sim_Wheat": taux_simul["Wheat"],
                    "Sim_Alfalfa": taux_simul["Alfalfa"],
                    
                    "Obs_Total_moy": row['Obs_Total_moy'], "Obs_Total_std": row['Obs_Total_std'],
                    "Obs_Center_moy": row['Obs_Center_moy'], "Obs_Center_std": row['Obs_Center_std'],
                    "Obs_Side_moy": row['Obs_Side_moy'], "Obs_Side_std": row['Obs_Side_std'],
                    "Obs_IR_moy": row['Obs_IR_moy'], "Obs_IR_std": row['Obs_IR_std'],
                })

    # ==========================================
    # 4. GRAPHIQUES D'ÉVOLUTION PAR ESPÈCE
    # ==========================================
    if resultats:
        df_res = pd.DataFrame(resultats)
        df_res.to_csv(os.path.join(DOSSIER_SAUVEGARDE, "Comparaison_Detaillee.csv"), index=False, sep=';', decimal=',')
        
        for simulation in df_res['Simulation'].unique():
            df_plot = df_res[df_res['Simulation'] == simulation].sort_values(by='DAS')
            is_wheatluz = "WheatLuz" in simulation
            
            fig, ax1 = plt.subplots(figsize=(11, 6))
            
            # --- TOTAL ---
            ax1.plot(df_plot['DAS'], df_plot['Sim_Total'], marker='o', linestyle='-', linewidth=2.5, color='black', label='Simulé - Total')
            ax1.errorbar(df_plot['DAS'], df_plot['Obs_Total_moy'], yerr=df_plot['Obs_Total_std'], fmt='s--', color='gray', label='Observé - Total')
            
            # --- BLÉ ---
            ax1.plot(df_plot['DAS'], df_plot['Sim_Wheat'], marker='^', linestyle='-', linewidth=2, color='green', label='Simulé - Blé')
            ax1.errorbar(df_plot['DAS'], df_plot['Obs_Center_moy'], yerr=df_plot['Obs_Center_std'], fmt='^--', color='limegreen', label='Observé - Center (Blé)')
            
            if not is_wheatluz:
                # Si culture pure de blé, on trace également le Side en tant que Blé
                ax1.errorbar(df_plot['DAS'], df_plot['Obs_Side_moy'], yerr=df_plot['Obs_Side_std'], fmt='x--', color='darkgreen', label='Observé - Side (Blé)')

            # --- LUZERNE ---
            if is_wheatluz:
                ax1.plot(df_plot['DAS'], df_plot['Sim_Alfalfa'], marker='v', linestyle='-', linewidth=2, color='purple', label='Simulé - Luzerne')
                ax1.errorbar(df_plot['DAS'], df_plot['Obs_Side_moy'], yerr=df_plot['Obs_Side_std'], fmt='v--', color='mediumorchid', label='Observé - Side (Luzerne)')
            
            ax1.set_title(f"Dynamique de Couverture par Zone/Espèce - {simulation}\n", fontweight='bold')
            ax1.set_xlabel("Jours après semis (DAS)", fontweight='bold')
            ax1.set_ylabel("Taux de couverture (0-1)", fontweight='bold')
            ax1.set_xticks(df_plot['DAS'].astype(int))
            
            ax1.legend(loc="upper left", bbox_to_anchor=(1.02, 1), borderaxespad=0.)
            ax1.grid(True, linestyle='--', alpha=0.5)
            
            ax2 = ax1.twiny()
            ax2.set_xlim(ax1.get_xlim())
            ax2.set_xticks(df_plot['DAS'].values)
            ax2.set_xticklabels(df_plot['Date'].values)
            ax2.set_xlabel("Date (JJ/MM)", fontweight='bold')

            plt.savefig(os.path.join(DOSSIER_SAUVEGARDE, f"Evolution_Detaillee_{simulation}.png"), dpi=150, bbox_inches='tight')
            plt.close()

        # ==========================================
        # 5. GRAPHIQUE RMSE
        # ==========================================
        df_rmse = df_res.copy()
        if DATE_LIMITE_RMSE:
            df_rmse = df_rmse[df_rmse['Date_dt'] <= pd.to_datetime(DATE_LIMITE_RMSE)]
            
        for simulation in df_rmse['Simulation'].unique():
            df_sim = df_rmse[df_rmse['Simulation'] == simulation]
            is_wheatluz = "WheatLuz" in simulation
            
            fig, ax = plt.subplots(figsize=(6, 6))
            
            # RMSE Total
            rmse_tot = np.sqrt(np.mean((df_sim['Sim_Total'] - df_sim['Obs_Total_moy'])**2))
            ax.scatter(df_sim['Obs_Total_moy'], df_sim['Sim_Total'], c='black', marker='o', label=f'Total (RMSE={rmse_tot:.3f})')
            
            # RMSE Blé (Simulé vs Center)
            rmse_ble = np.sqrt(np.mean((df_sim['Sim_Wheat'] - df_sim['Obs_Center_moy'])**2))
            ax.scatter(df_sim['Obs_Center_moy'], df_sim['Sim_Wheat'], c='green', marker='^', label=f'Blé vs Center (RMSE={rmse_ble:.3f})')
            
            # RMSE Luzerne (Simulé vs Side)
            if is_wheatluz:
                rmse_luz = np.sqrt(np.mean((df_sim['Sim_Alfalfa'] - df_sim['Obs_Side_moy'])**2))
                ax.scatter(df_sim['Obs_Side_moy'], df_sim['Sim_Alfalfa'], c='purple', marker='v', label=f'Luzerne vs Side (RMSE={rmse_luz:.3f})')

            ax.plot([0, 1], [0, 1], 'r--', zorder=1)
            ax.set_title(f"RMSE par zone - {simulation}", fontweight='bold')
            ax.set_xlabel("Observé (CSV)")
            ax.set_ylabel("Simulé (VTK)")
            ax.set_xlim(0, max(df_sim['Obs_Total_moy'].max(), df_sim['Sim_Total'].max()) * 1.1)
            ax.set_ylim(0, max(df_sim['Obs_Total_moy'].max(), df_sim['Sim_Total'].max()) * 1.1)
            ax.legend(loc='lower right')
            ax.grid(True, linestyle=':', alpha=0.6)
            
            plt.savefig(os.path.join(DOSSIER_SAUVEGARDE, f"RMSE_Detaille_{simulation}.png"), dpi=150, bbox_inches='tight')
            plt.close()

    print("\n✅ Traitement terminé. Graphiques enregistrés dans le dossier spécifié.")