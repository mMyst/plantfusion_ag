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
    "WheatLuz_LUBBAC_emergence_adel": [2,3,5,9], #[2, 3, 4, 5, 6, 7, 8, 9], #
    "Wheat_LUBBAC_emergence_adel": [1],
}

# Mise à jour avec le nouveau fichier
FICHIER_CSV = "resultats_taux_couv_detaille.csv"

DATE_SEMIS = pd.to_datetime("2024-11-18") # Jour 0 (DAS 0)
DATE_LIMITE_RMSE = None #"2025-03-15"  

X_DOMAINE = 0.3      
RESOLUTION = 0.001   
HAUTEUR_CIBLE = 600  
LARGEUR_CIBLE = 1150 
DECALAGE_X = 50      

DOSSIER_SAUVEGARDE = r"C:\Users\agrumel\Documents\Data\Simul LUBBAC\Sauvergardes_Figures\latest\Taux_couverture_par_zone"

# ==========================================
# 2. CALCUL DU TAUX DE COUVERTURE DEPUIS LE VTK
# ==========================================
def generer_vtk_taux_couverture_spatial(dossier_vtk, doy, X_dom, res, h_cible, w_cible):
    """
    Lit le fichier VTK, rastérise la scène globale puis calcule la couverture de
    chaque zone PAR RAPPORT AU CARRÉ ENTIER (l'image complète).
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
    except Exception:
        return None

    taille_base = int(X_dom / res)
    z_buffer_base = np.full((taille_base, taille_base), -np.inf)
    
    def f_bary(v0, v1, p):
        return (v1[0] - v0[0]) * (p[1] - v0[1]) - (v1[1] - v0[1]) * (p[0] - v0[0])

    for dx in [-X_dom, 0.0, X_dom]:
        for dy in [-X_dom, 0.0, X_dom]:
            pts_dec = pts.copy()
            pts_dec[:, 0] += dx
            pts_dec[:, 1] += dy
            
            for tri in cells:
                v0, v1, v2 = pts_dec[tri[0]], pts_dec[tri[1]], pts_dec[tri[2]]
                v0_p, v1_p, v2_p = v0[:2]/res, v1[:2]/res, v2[:2]/res
                
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

    z_buffer_base = np.where(z_buffer_base == -np.inf, 0.0, 1.0)
    z_buffer_base = np.roll(z_buffer_base, shift=DECALAGE_X, axis=1)
    
    reps_y = int(np.ceil(h_cible / taille_base))
    reps_x = int(np.ceil(w_cible / taille_base))
    
    img_tuilee = np.tile(z_buffer_base, (reps_y, reps_x))
    img_tuilee = img_tuilee[:h_cible, :w_cible] 
    
    # --- MASQUAGE SPATIAL ---
    cols = np.arange(w_cible)
    rem = cols % taille_base
    
    col_mask_side = (rem < 50) | (rem >= 250)
    col_mask_center = (rem >= 100) & (rem < 200)
    col_mask_ir = ((rem >= 50) & (rem < 100)) | ((rem >= 200) & (rem < 250))
    
    mask_side = np.tile(col_mask_side, (h_cible, 1))
    mask_center = np.tile(col_mask_center, (h_cible, 1))
    mask_ir = np.tile(col_mask_ir, (h_cible, 1))
    
    # CALCUL PAR RAPPORT AU CARRE ENTIER (Nombre de pixels total de l'image)
    total_pixels = img_tuilee.size
    
    return {
        "Total": np.sum(img_tuilee) / total_pixels,
        "Side": np.sum(img_tuilee[mask_side]) / total_pixels,
        "Center": np.sum(img_tuilee[mask_center]) / total_pixels,
        "IR": np.sum(img_tuilee[mask_ir]) / total_pixels
    }


# ==========================================
# 3. EXTRACTION DU CSV ET COMPARAISON
# ==========================================
if __name__ == "__main__":
    print("=======================================")
    print("COMPARAISON DÉTAILLÉE PAR ZONE (Ratio / Carré entier)")
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
        
        if not os.path.exists(dossier_vtk):
            continue
            
        print(f"\n>> Traitement : {simulation} (Bacs: {bacs_reels})")
        df_simul = df_csv[df_csv['Bac'].isin(bacs_reels)]
        
        if df_simul.empty: continue
            
        # Agrégation en utilisant les colonnes avec "_Carre"
        df_obs_stats = df_simul.groupby(['DAS', 'DOY', 'Date_dt', 'Date_Affiche']).agg(
            Obs_Total_moy=('Moyenne_Total', 'mean'), Obs_Total_std=('SD_Total', 'mean'),
            Obs_Side_moy=('Moyenne_Side_Carre', 'mean'), Obs_Side_std=('SD_Side_Carre', 'mean'),
            Obs_Center_moy=('Moyenne_Center_Carre', 'mean'), Obs_Center_std=('SD_Center_Carre', 'mean'),
            Obs_IR_moy=('Moyenne_IR_Carre', 'mean'), Obs_IR_std=('SD_IR_Carre', 'mean')
        ).reset_index().sort_values(by='DAS')

        for _, row in df_obs_stats.iterrows():
            taux_simul = generer_vtk_taux_couverture_spatial(dossier_vtk, int(row['DOY']), X_DOMAINE, RESOLUTION, HAUTEUR_CIBLE, LARGEUR_CIBLE)
            
            if taux_simul is not None:
                diff_tot = taux_simul["Total"] - row['Obs_Total_moy']
                print(f"  - DAS {int(row['DAS']):3d} | Sim Total: {taux_simul['Total']:.3f} | Obs Total: {row['Obs_Total_moy']:.3f} | Diff: {diff_tot:+.3f}")
                
                resultats.append({
                    "Simulation": simulation,
                    "DAS": int(row['DAS']),
                    "Date_dt": row['Date_dt'],
                    "Date": row['Date_Affiche'],
                    
                    "Sim_Total": taux_simul["Total"],
                    "Sim_Center": taux_simul["Center"],
                    "Sim_Side": taux_simul["Side"],
                    "Sim_IR": taux_simul["IR"],
                    
                    "Obs_Total_moy": row['Obs_Total_moy'], "Obs_Total_std": row['Obs_Total_std'],
                    "Obs_Center_moy": row['Obs_Center_moy'], "Obs_Center_std": row['Obs_Center_std'],
                    "Obs_Side_moy": row['Obs_Side_moy'], "Obs_Side_std": row['Obs_Side_std'],
                    "Obs_IR_moy": row['Obs_IR_moy'], "Obs_IR_std": row['Obs_IR_std'],
                })

    # ==========================================
    # 4. GRAPHIQUES D'ÉVOLUTION PAR ZONE SPATIALE
    # ==========================================
    if resultats:
        df_res = pd.DataFrame(resultats)
        df_res.to_csv(os.path.join(DOSSIER_SAUVEGARDE, "Comparaison_Spatiale_Zones_Carre.csv"), index=False, sep=';', decimal=',')
        
        for simulation in df_res['Simulation'].unique():
            df_plot = df_res[df_res['Simulation'] == simulation].sort_values(by='DAS')
            is_wheatluz = "WheatLuz" in simulation
            
            lbl_center = "Center (Blé)"
            lbl_side = "Side (Luzerne)" if is_wheatluz else "Side (Blé)"
            
            fig, ax1 = plt.subplots(figsize=(12, 7))
            
            # --- TOTAL ---
            ax1.plot(df_plot['DAS'], df_plot['Sim_Total'], marker='o', linestyle='-', linewidth=2.5, color='black', label='Simulé - Total')
            ax1.errorbar(df_plot['DAS'], df_plot['Obs_Total_moy'], yerr=df_plot['Obs_Total_std'], fmt='s--', color='gray', label='Observé - Total')
            
            # --- CENTER ---
            ax1.plot(df_plot['DAS'], df_plot['Sim_Center'], marker='^', linestyle='-', linewidth=2, color='green', label=f'Simulé - {lbl_center}')
            ax1.errorbar(df_plot['DAS'], df_plot['Obs_Center_moy'], yerr=df_plot['Obs_Center_std'], fmt='^--', color='limegreen', label=f'Observé - {lbl_center}')
            
            # --- SIDE ---
            color_sim_side = 'purple' if is_wheatluz else 'darkgreen'
            color_obs_side = 'mediumorchid' if is_wheatluz else 'forestgreen'
            
            ax1.plot(df_plot['DAS'], df_plot['Sim_Side'], marker='v', linestyle='-', linewidth=2, color=color_sim_side, label=f'Simulé - {lbl_side}')
            ax1.errorbar(df_plot['DAS'], df_plot['Obs_Side_moy'], yerr=df_plot['Obs_Side_std'], fmt='v--', color=color_obs_side, label=f'Observé - {lbl_side}')
                
            # --- INTER-RANG (IR) ---
            ax1.plot(df_plot['DAS'], df_plot['Sim_IR'], marker='d', linestyle='-', linewidth=2, color='saddlebrown', label='Simulé - Inter-rang (IR)')
            ax1.errorbar(df_plot['DAS'], df_plot['Obs_IR_moy'], yerr=df_plot['Obs_IR_std'], fmt='d--', color='peru', label='Observé - Inter-rang (IR)')
            
            ax1.set_title(f"Dynamique de Couverture par Zone Spatiale (Rapportée au carré) - {simulation}\n", fontweight='bold')
            ax1.set_xlabel("Jours après semis (DAS)", fontweight='bold')
            ax1.set_ylabel("Taux de couverture (0-1)", fontweight='bold')
            ax1.set_xticks(df_plot['DAS'].astype(int))
            
            ax1.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), borderaxespad=0.)
            ax1.grid(True, linestyle='--', alpha=0.5)
            
            ax2 = ax1.twiny()
            ax2.set_xlim(ax1.get_xlim())
            ax2.set_xticks(df_plot['DAS'].values)
            ax2.set_xticklabels(df_plot['Date'].values)
            ax2.set_xlabel("Date (JJ/MM)", fontweight='bold')

            plt.savefig(os.path.join(DOSSIER_SAUVEGARDE, f"Evolution_Zones_Carre_{simulation}.png"), dpi=150, bbox_inches='tight')
            plt.close()

        # ==========================================
        # 5. GRAPHIQUE RMSE PAR ZONE SPATIALE
        # ==========================================
        df_rmse = df_res.copy()
        if DATE_LIMITE_RMSE:
            df_rmse = df_rmse[df_rmse['Date_dt'] <= pd.to_datetime(DATE_LIMITE_RMSE)]
            
        for simulation in df_rmse['Simulation'].unique():
            df_sim = df_rmse[df_rmse['Simulation'] == simulation]
            is_wheatluz = "WheatLuz" in simulation
            
            lbl_center = "Center (Blé)"
            lbl_side = "Side (Luzerne)" if is_wheatluz else "Side (Blé)"
            color_side = 'purple' if is_wheatluz else 'darkgreen'
            
            fig, ax = plt.subplots(figsize=(7, 7))
            
            # Total
            rmse_tot = np.sqrt(np.mean((df_sim['Sim_Total'] - df_sim['Obs_Total_moy'])**2))
            ax.scatter(df_sim['Obs_Total_moy'], df_sim['Sim_Total'], c='black', marker='o', label=f'Total (RMSE={rmse_tot:.3f})')
            
            # Center
            rmse_center = np.sqrt(np.mean((df_sim['Sim_Center'] - df_sim['Obs_Center_moy'])**2))
            ax.scatter(df_sim['Obs_Center_moy'], df_sim['Sim_Center'], c='green', marker='^', label=f'{lbl_center} (RMSE={rmse_center:.3f})')
            
            # Side
            rmse_side = np.sqrt(np.mean((df_sim['Sim_Side'] - df_sim['Obs_Side_moy'])**2))
            ax.scatter(df_sim['Obs_Side_moy'], df_sim['Sim_Side'], c=color_side, marker='v', label=f'{lbl_side} (RMSE={rmse_side:.3f})')

            # IR
            rmse_ir = np.sqrt(np.mean((df_sim['Sim_IR'] - df_sim['Obs_IR_moy'])**2))
            ax.scatter(df_sim['Obs_IR_moy'], df_sim['Sim_IR'], c='saddlebrown', marker='d', label=f'Inter-rang (RMSE={rmse_ir:.3f})')

            ax.plot([0, 1], [0, 1], 'r--', zorder=1)
            ax.set_title(f"RMSE par zone (Rapportée au carré) - {simulation}", fontweight='bold')
            ax.set_xlabel("Observé (CSV)")
            ax.set_ylabel("Simulé (VTK)")
            
            max_val = max(df_sim[['Obs_Total_moy', 'Obs_Center_moy', 'Obs_Side_moy', 'Obs_IR_moy']].max().max(), 
                          df_sim[['Sim_Total', 'Sim_Center', 'Sim_Side', 'Sim_IR']].max().max())
            ax.set_xlim(0, max_val * 1.1)
            ax.set_ylim(0, max_val * 1.1)
            ax.set_aspect('equal', adjustable='box')
            ax.legend(loc='lower right', fontsize='small')
            ax.grid(True, linestyle=':', alpha=0.6)
            
            plt.savefig(os.path.join(DOSSIER_SAUVEGARDE, f"RMSE_Zones_Carre_{simulation}.png"), dpi=150, bbox_inches='tight')
            plt.close()

    print("\n✅ Traitement spatial terminé. Graphiques enregistrés dans le dossier spécifié.")