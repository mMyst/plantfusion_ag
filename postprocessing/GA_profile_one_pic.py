import os
import re
import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize

def calculer_profil_z(sous_maillage, dz):
    if sous_maillage.n_cells == 0:
        return np.array([]), np.array([])
    centres = sous_maillage.cell_centers().points
    z_centres = centres[:, 2]
    surfaces = sous_maillage.cell_data['Area']
    z_min, z_max = z_centres.min(), z_centres.max()
    bins = np.arange(z_min, z_max + dz, dz)
    surface_par_z, _ = np.histogram(z_centres, bins=bins, weights=surfaces)
    z_moyens = bins[:-1] + (dz / 2.0)
    return z_moyens, surface_par_z

def process_vtk_canopy(file_path, dz):
    mesh = pv.read(file_path)
    mesh_ble = mesh.threshold([0.9, 1.1], scalars='VegetationType')
    mesh_luzerne = mesh.threshold([-0.1, 0.1], scalars='VegetationType')
    
    if mesh_ble.n_cells > 0:
        mesh_ble = mesh_ble.compute_cell_sizes(length=False, volume=False)
    if mesh_luzerne.n_cells > 0:
        mesh_luzerne = mesh_luzerne.compute_cell_sizes(length=False, volume=False)
        
    z_ble, surf_ble = calculer_profil_z(mesh_ble, dz)
    z_luzerne, surf_luzerne = calculer_profil_z(mesh_luzerne, dz)
    return z_ble, surf_ble, z_luzerne, surf_luzerne

def main():
    # ================= CONFIGURATION =================
    doy_start = 14                    # Premier DOY à traiter
    doy_end = 125                     # Dernier DOY à traiter (inclus)
    vtk_dir = r'C:\Users\agrumel\Code\Python_Ecophy\plantfusion_ag\outputs\WheatLuzREG_LUBBAC_emergence_na\light\vtk'                # Dossier contenant les VTK
    out_dir = "./profils_canopee/all"   # Dossier où les graphiques seront enregistrés
    dz = 0.05                         # Incrément de hauteur z (en mètres)                    # Incrément de hauteur z (en mètres)
    # =================================================


    os.makedirs(out_dir, exist_ok=True)
    
    doys_list = range(doy_start, doy_end + 1)
    print(f"Recherche des fichiers VTK pour les DOY de {doy_start} à {doy_end}...")
    
    resultats_ble = {}
    resultats_luzerne = {}
    
    try:
        tous_les_fichiers = os.listdir(vtk_dir)
    except FileNotFoundError:
        print(f"Erreur : Le dossier {vtk_dir} est introuvable.")
        return
    
    fichiers_traites_count = 0

    for doy in doys_list:
        doy_str = str(doy).zfill(3)
        # Utilisation de l'expression régulière demandée
        motif = re.compile(rf"^.*_\d+{doy_str}13\.vtk$")
        
        fichiers_match = [f for f in tous_les_fichiers if motif.match(f)]
        
        if not fichiers_match:
            continue
            
        fichier_cible = fichiers_match[0]
        chemin_complet = os.path.join(vtk_dir, fichier_cible)
        
        print(f" [*] Traitement du DOY {doy_str} -> {fichier_cible}")
        z_ble, surf_ble, z_luzerne, surf_luzerne = process_vtk_canopy(chemin_complet, dz)
        
        # Stocker les valeurs pour le graphique de synthèse
        if len(z_ble) > 0:
            resultats_ble[doy] = (z_ble, surf_ble)
        if len(z_luzerne) > 0:
            resultats_luzerne[doy] = (z_luzerne, surf_luzerne)
            
        fichiers_traites_count += 1
            
    if fichiers_traites_count == 0:
        print("Aucun fichier VTK correspondant n'a été trouvé dans la plage spécifiée.")
        return
        
    print(f"{fichiers_traites_count} fichiers traités. Génération du graphique de synthèse...")
    fig, axs = plt.subplots(2, 2, figsize=(16, 12), sharey=True)
    
    doys_traites = sorted(list(set(list(resultats_ble.keys()) + list(resultats_luzerne.keys()))))
    norm = Normalize(vmin=min(doys_traites), vmax=max(doys_traites))
    cmap = cm.viridis
    
    for d in doys_traites:
        color = cmap(norm(d))
        
        if d in resultats_ble:
            zb, sb = resultats_ble[d]
            axs[0, 0].plot(sb, zb, color=color, alpha=0.8, linewidth=1.5)
            axs[1, 0].plot(np.cumsum(sb), zb, color=color, alpha=0.8, linewidth=1.5)
            
        if d in resultats_luzerne:
            zl, sl = resultats_luzerne[d]
            axs[0, 1].plot(sl, zl, color=color, alpha=0.8, linewidth=1.5)
            axs[1, 1].plot(np.cumsum(sl), zl, color=color, alpha=0.8, linewidth=1.5)
            
    axs[0, 0].set_title("Blé - Densité de surface")
    axs[0, 1].set_title("Luzerne - Densité de surface")
    axs[1, 0].set_title("Blé - Surface Cumulée")
    axs[1, 1].set_title("Luzerne - Surface Cumulée")
    
    for ax in axs.flat:
        ax.set_ylabel("Hauteur Z (m)")
        ax.grid(True, linestyle='--', alpha=0.5)
        
    axs[1, 0].set_xlabel("Surface cumulée ($m^2$)")
    axs[1, 1].set_xlabel("Surface cumulée ($m^2$)")
    axs[0, 0].set_xlabel("Surface ($m^2$)")
    axs[0, 1].set_xlabel("Surface ($m^2$)")
    
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axs, orientation='vertical', fraction=0.02, pad=0.04)
    cbar.set_label('DOY (Day Of Year)')
    
    plt.savefig(os.path.join(out_dir, "Graphique_Synthese.png"), dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"Terminé avec succès. Graphique de synthèse sauvegardé dans le dossier : {out_dir}")

if __name__ == '__main__':
    main()

