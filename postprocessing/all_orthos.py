import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from PIL import Image
import warnings

# ==========================================
# 1. PARAMÈTRES ET CONFIGURATION
# ==========================================
CHEMIN_BASE = "C:/Users/agrumel/Documents/Data/Simul LUBBAC"
CHEMIN_SUFFIXE = "light/vtk"

CONFIG_ROUTAGE = {
    "WheatLuz_LUBBAC_emergence_adel": [2, 5],
    "Wheat_LUBBAC_emergence_adel": [1],
}

FICHIER_EXCEL = "//pnas1.stockage.inra.fr/nap-urp3f/Projets/2024_Lubbac_GLouarn/Orthophoto/correspondances_date_DOY.xlsx"

X_DOMAINE = 0.3      
RESOLUTION = 0.001   
HAUTEUR_CIBLE = 600  
LARGEUR_CIBLE = 1150 
DECALAGE_X = 50      

DOSSIER_SAUVEGARDE = r"C:\Users\agrumel\Documents\Data\Simul LUBBAC\Sauvergardes_Figures\latest\Resultats_Comparaisons"

# ==========================================
# 2. GÉNÉRATION DE LA TUILE VTK
# ==========================================
def generer_vtk_tuile(dossier_vtk, f_cible, doy, X_dom, res, h_cible, w_cible):
    if not os.path.exists(dossier_vtk):
        return np.full((h_cible, w_cible), np.nan)

    fichiers = os.listdir(dossier_vtk)
    chemin_vtk = None
    doy_str = str(doy).zfill(3) 
    
    motif = rf"^.*_\d+{doy_str}13\.vtk$" 
    for f in fichiers:
        if re.match(motif, f):
            chemin_vtk = os.path.join(dossier_vtk, f)
            break
            
    if not chemin_vtk:
        return np.full((h_cible, w_cible), np.nan)
        
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
        return np.full((h_cible, w_cible), np.nan)

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

    if f_cible == "taux_couv":
        z_buffer_base = np.where(z_buffer_base == -np.inf, 0.0, 1.0)
    else:
        z_buffer_base[z_buffer_base == -np.inf] = np.nan

    z_buffer_base = np.roll(z_buffer_base, shift=DECALAGE_X, axis=1)
    
    reps_y = int(np.ceil(h_cible / taille_base))
    reps_x = int(np.ceil(w_cible / taille_base))
    
    img_tuilee = np.tile(z_buffer_base, (reps_y, reps_x))
    img_tuilee = img_tuilee[:h_cible, :w_cible] 
    
    return img_tuilee

# ==========================================
# 3. EXTRACTION DES RÉELLES & UTILITAIRES
# ==========================================
def traiter_images_reelles(doy, f_cible, bacs_cible):
    try:
        df = pd.read_excel(FICHIER_EXCEL)
    except Exception:
        return None, None

    df_doy = df[df['DOY'].astype(str) == str(doy)]
    df_principal = df_doy[df_doy['type_fichier'] == f_cible]
    df_masque = df_doy[df_doy['type_fichier'] == "taux_couv"] if f_cible == "hauteur" else pd.DataFrame()

    images_chargees = []
    
    for bac in bacs_cible:
        colonne = f"selec_veg_bac{bac}"
        
        if colonne in df_principal.columns:
            chemins = df_principal[colonne].dropna().values
            for chemin in chemins:
                if isinstance(chemin, str) and chemin.strip() != "":
                    chemin = chemin.replace('\\', '/')
                    if os.path.exists(chemin):
                        img = np.array(Image.open(chemin)).astype(float)
                        if img.ndim == 3: 
                            img = img[:, :, 0]
                            
                        if f_cible == "hauteur":
                            img = np.where(img <= 10000, np.nan, (img - 10000.0) / 1000.0)
                            if not df_masque.empty and colonne in df_masque.columns:
                                chemins_masque = df_masque[colonne].dropna().values
                                if len(chemins_masque) > 0:
                                    chemin_masque = chemins_masque[0].replace('\\', '/')
                                    if os.path.exists(chemin_masque):
                                        img_m = np.array(Image.open(chemin_masque)).astype(float)
                                        if img_m.ndim == 3: img_m = img_m[:, :, 0]
                                        img = np.where(img_m > 0, img, np.nan)

                        elif f_cible == "taux_couv":
                            img = np.where(img > 0, 1.0, 0.0)
                            
                        images_chargees.append(img)
    
    if not images_chargees:
        return np.full((HAUTEUR_CIBLE, LARGEUR_CIBLE), np.nan), np.full((HAUTEUR_CIBLE, LARGEUR_CIBLE), np.nan)

    pile_images = np.stack(images_chargees, axis=0)
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        img_moyenne = np.nanmean(pile_images, axis=0)
        img_sd = np.nanstd(pile_images, axis=0)
    
    return img_moyenne, img_sd

def calcul_moyenne_top_1_pourcent(img):
    """Calcule la moyenne des 1% des valeurs les plus hautes pour chaque colonne X."""
    W = img.shape[1]
    profil = np.full(W, np.nan)
    for x in range(W):
        colonne = img[:, x]
        valides = colonne[~np.isnan(colonne)] # On exclut les NaN
        if len(valides) > 0:
            # On prend au moins 1 pixel, sinon les 1% des pixels valides
            k = max(1, int(np.ceil(0.01 * len(valides))))
            profil[x] = np.mean(np.sort(valides)[-k:])
    return profil

# ==========================================
# 4. ROUTINE DE COMPARAISON & SAUVEGARDE
# ==========================================
def routine_comparaison(nom_simul, f_cible, doy_cible):
    bacs_reels = CONFIG_ROUTAGE.get(nom_simul, [])
    if not bacs_reels:
        return False

    dossier_vtk = os.path.join(CHEMIN_BASE, nom_simul, CHEMIN_SUFFIXE).replace('\\', '/')
    
    img_vtk = generer_vtk_tuile(dossier_vtk, f_cible, doy_cible, X_DOMAINE, RESOLUTION, HAUTEUR_CIBLE, LARGEUR_CIBLE)
    
    if np.all(np.isnan(img_vtk)):
        return False

    img_reelle_moy, img_reelle_sd = traiter_images_reelles(doy_cible, f_cible, bacs_reels)
    
    if np.all(np.isnan(img_reelle_moy)):
        print(f"  -> Images réelles introuvables pour le DOY {doy_cible}, graphe ignoré.")
        return False
    
    cmap_main = 'viridis' if f_cible == "hauteur" else 'Greens'
    label_cb = 'Hauteur Z (m)' if f_cible == "hauteur" else 'Taux de couverture (0-1)'
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        profil_vtk_moy = np.nanmean(img_vtk, axis=0)
        profil_reelle_moy = np.nanmean(img_reelle_moy, axis=0)
        profil_vtk_max = np.nanmax(img_vtk, axis=0)
        profil_reelle_max = np.nanmax(img_reelle_moy, axis=0)
        
        if f_cible == "hauteur":
            profil_vtk_top1 = calcul_moyenne_top_1_pourcent(img_vtk)
            profil_reelle_top1 = calcul_moyenne_top_1_pourcent(img_reelle_moy)

    # --- AJUSTEMENT DYNAMIQUE DES GRAPHIQUES ---
    nb_axes = 6 if f_cible == "hauteur" else 5
    hauteur_fig = 12 if f_cible == "hauteur" else 10
    ratios = [1.5, 1.5, 1.5, 1, 1, 1] if f_cible == "hauteur" else [1.5, 1.5, 1.5, 1, 1]

    fig, axes = plt.subplots(nb_axes, 1, figsize=(6, hauteur_fig), constrained_layout=True, sharex=True, 
                             gridspec_kw={'height_ratios': ratios})
    
    im0 = axes[0].imshow(img_vtk, cmap=cmap_main, origin='lower', aspect='auto')
    axes[0].set_title(f"VTK Tuilé (Simulation : {nom_simul} - DOY : {doy_cible})", fontweight='bold', fontsize=10)
    fig.colorbar(im0, ax=axes[0], label=label_cb, pad=0.01)
    
    im1 = axes[1].imshow(img_reelle_moy, cmap=cmap_main, origin='lower', aspect='auto')
    axes[1].set_title(f"Moyenne Réelles (Bacs {bacs_reels})", fontweight='bold', fontsize=10)
    fig.colorbar(im1, ax=axes[1], label=label_cb, pad=0.01)
    
    im2 = axes[2].imshow(img_reelle_sd, cmap='plasma', origin='lower', aspect='auto')
    axes[2].set_title("Écart-Type (SD) inter-bacs", fontweight='bold', fontsize=10)
    fig.colorbar(im2, ax=axes[2], label='Variation (SD)', pad=0.01)
    
    for ax in axes[:3]:
        ax.set_ylabel("Pixels Y (mm)", fontsize=10)
        ax.set_facecolor('lightgrey' if f_cible == "hauteur" else 'white')
        
    x_valeurs = np.arange(LARGEUR_CIBLE)
    
    # 3. MOYENNE
    axes[3].plot(x_valeurs, profil_vtk_moy, label='VTK Tuilé', color='blue', linewidth=2)
    axes[3].plot(x_valeurs, profil_reelle_moy, label='Réelles', color='orange', linewidth=2, alpha=0.8)
    axes[3].set_title("Comparaison des profils MOYENS", fontweight='bold')
    axes[3].set_ylabel(f"Moyenne {label_cb}", fontsize=6)
    axes[3].legend(loc='upper right')
    axes[3].grid(True, linestyle='--', alpha=0.7) 
    
    # 4. MAXIMUM
    axes[4].plot(x_valeurs, profil_vtk_max, label='VTK Tuilé', color='blue', linewidth=2)
    axes[4].plot(x_valeurs, profil_reelle_max, label='Réelles', color='orange', linewidth=2, alpha=0.8)
    axes[4].set_title("Comparaison des profils MAXIMUM", fontweight='bold')
    axes[4].set_ylabel(f"Max {label_cb}", fontsize=6)
    axes[4].legend(loc='upper right')
    axes[4].grid(True, linestyle='--', alpha=0.7) 
    
    # 5. TOP 1% (Si hauteur)
    if f_cible == "hauteur":
        axes[5].plot(x_valeurs, profil_vtk_top1, label='VTK Tuilé', color='blue', linewidth=2)
        axes[5].plot(x_valeurs, profil_reelle_top1, label='Réelles', color='orange', linewidth=2, alpha=0.8)
        axes[5].set_title("Moyenne des 1% les plus hauts (Enveloppe canopée lissée)", fontweight='bold')
        axes[5].set_xlabel("Pixels X (mm)", fontsize=10) 
        axes[5].set_ylabel(f"Top 1% {label_cb}", fontsize=6)
        axes[5].legend(loc='upper right')
        axes[5].grid(True, linestyle='--', alpha=0.7) 
    else:
        axes[4].set_xlabel("Pixels X (mm)", fontsize=10) 

    # --- SAUVEGARDE ENREGISTREMENT ---
    nom_fichier = f"{nom_simul}_{f_cible}_DOY{doy_cible}.png"
    chemin_complet = os.path.join(DOSSIER_SAUVEGARDE, nom_fichier)
    
    plt.savefig(chemin_complet, dpi=150)
    plt.close(fig) 
    
    print(f"✅ Graphe sauvegardé : {nom_fichier}")
    return True

# ==========================================
# 5. SCRIPT D'AUTOMATISATION DE MASSE
# ==========================================
if __name__ == "__main__":
    print("=======================================")
    print("LANCEMENT DE LA GÉNÉRATION DE MASSE")
    print("=======================================")
    
    if not os.path.exists(DOSSIER_SAUVEGARDE):
        os.makedirs(DOSSIER_SAUVEGARDE)
        print(f"Création du dossier de sauvegarde : {DOSSIER_SAUVEGARDE}/")

    modes_analyse = ["taux_couv", "hauteur"]

    try:
        df_excel = pd.read_excel(FICHIER_EXCEL)
        doys_excel = sorted([str(int(d)) for d in df_excel['DOY'].dropna().unique()])
        print(f"✅ {len(doys_excel)} DOYs uniques identifiés dans l'Excel : {', '.join(doys_excel)}")
    except Exception as e:
        print(f"🛑 Erreur fatale : Impossible de lire l'Excel pour extraire les DOYs ({e})")
        exit()

    for simulation in CONFIG_ROUTAGE.keys():
        dossier_vtk = os.path.join(CHEMIN_BASE, simulation, CHEMIN_SUFFIXE).replace('\\', '/')
        
        if not os.path.exists(dossier_vtk):
            print(f"\n⚠️ IGNORÉ : Le dossier pour la simulation '{simulation}' n'existe pas : {dossier_vtk}")
            continue
            
        print(f"\n>> Traitement de la simulation : {simulation}")

        for mode in modes_analyse:
            print(f"\n--- Traitement Mode : {mode.upper()} ---")
            for doy in doys_excel:
                routine_comparaison(simulation, mode, doy)

    print("\n=======================================")
    print("FIN DU TRAITEMENT BATCH !")
    print("=======================================")