import os
import glob
import csv
import re
import numpy as np
from PIL import Image

# ================= CONFIGURATION =================
DOSSIER_IMAGES = r"\\147.100.80.181\Phenotypage\Data\Lubac_2024\Qgis\Image\mask"
FICHIER_CSV = "resultats_taux_couv_detaille.csv"
# Formats d'images autorisés (ignore automatiquement les .csv et autres fichiers)
EXTENSIONS_IMAGES = (".jpg", ".jpeg", ".png", ".tif", ".tiff", ".bmp")
# =================================================
CLES_METRIQUES = [
    "total",
    "side_local", "ir_local", "center_local",
    "side_carre", "ir_carre", "center_carre",
    "part_veg_side", "part_veg_ir", "part_veg_center"
]

def analyser_carre(crop_blanc, crop_valide):
    """
    crop_blanc  : matrice booléenne des pixels blancs non-transparents
    crop_valide : matrice booléenne des pixels opaques (sol ou végétation)
    """
    nb_valide_carre = np.count_nonzero(crop_valide)
    nb_blanc_carre = np.count_nonzero(crop_blanc)

    # Découpage des masques selon les sous-zones x (0-300)
    # Side : 0-50 et 250-300
    side_b = np.concatenate([crop_blanc[:, :50], crop_blanc[:, 250:300]], axis=1)
    side_v = np.concatenate([crop_valide[:, :50], crop_valide[:, 250:300]], axis=1)

    # IR : 50-100 et 200-250
    ir_b = np.concatenate([crop_blanc[:, 50:100], crop_blanc[:, 200:250]], axis=1)
    ir_v = np.concatenate([crop_valide[:, 50:100], crop_valide[:, 200:250]], axis=1)

    # Center : 100-200
    center_b = crop_blanc[:, 100:200]
    center_v = crop_valide[:, 100:200]

    # Comptages
    nb_side_b, nb_side_v = np.count_nonzero(side_b), np.count_nonzero(side_v)
    nb_ir_b, nb_ir_v = np.count_nonzero(ir_b), np.count_nonzero(ir_v)
    nb_center_b, nb_center_v = np.count_nonzero(center_b), np.count_nonzero(center_v)

    def ratio(numerateur, denominateur):
        return (numerateur / denominateur) if denominateur > 0 else 0.0

    return {
        # 1. Taux global du carré (sur la surface opaque réelle)
        "total": ratio(nb_blanc_carre, nb_valide_carre),

        # 2. Taux locaux par zone (par rapport aux pixels opaques de chaque zone)
        "side_local": ratio(nb_side_b, nb_side_v),
        "ir_local": ratio(nb_ir_b, nb_ir_v),
        "center_local": ratio(nb_center_b, nb_center_v),

        # 3. Taux rapportés à la surface opaque totale du carré
        "side_carre": ratio(nb_side_b, nb_valide_carre),
        "ir_carre": ratio(nb_ir_b, nb_valide_carre),
        "center_carre": ratio(nb_center_b, nb_valide_carre),

        # 4. Part de la végétation présente dans chaque zone
        "part_veg_side": ratio(nb_side_b, nb_blanc_carre),
        "part_veg_ir": ratio(nb_ir_b, nb_blanc_carre),
        "part_veg_center": ratio(nb_center_b, nb_blanc_carre),
    }

def generer_en_tetes():
    en_tetes = ["Nom_Photo", "Date", "Bac"]
    for num in (1, 2, 3):
        en_tetes.extend([
            f"C{num}_Total",
            f"C{num}_Side_Local", f"C{num}_IR_Local", f"C{num}_Center_Local",
            f"C{num}_Side_sur_Carre", f"C{num}_IR_sur_Carre", f"C{num}_Center_sur_Carre",
            f"C{num}_Part_Veg_Side", f"C{num}_Part_Veg_IR", f"C{num}_Part_Veg_Center"
        ])
    for cle in CLES_METRIQUES:
        nom_base = cle.replace("_", " ").title().replace(" ", "_")
        en_tetes.extend([f"Moyenne_{nom_base}", f"SD_{nom_base}"])
    return en_tetes

def traiter_images_dossier():
    # Sélection exclusive des fichiers PNG
    motif = os.path.join(DOSSIER_IMAGES, "Lubac_*.png")
    fichiers = glob.glob(motif)
    en_tetes = generer_en_tetes()
    donnees_csv = []

    # Regex acceptant indifféremment '_veg' ou '_mask' à la fin du nom
    pattern_nom = re.compile(r"(\d{4})_(\d{2})_(\d{2})_bac_([^_]+)_(?:veg|mask)$")

    for chemin in sorted(fichiers):
        nom_fichier = os.path.basename(chemin)
        nom_sans_ext, ext = os.path.splitext(nom_fichier)

        if not (nom_sans_ext.endswith("_veg") or nom_sans_ext.endswith("_mask")):
            continue

        match = pattern_nom.search(nom_sans_ext)
        if match:
            annee, mois, jour, bac = match.groups()
            date_formatee = f"{jour}/{mois}/{annee}"
        else:
            date_formatee = ""
            bac = ""

        try:
            with Image.open(chemin) as img:
                # Normalisation en RGBA pour isoler proprement le canal Alpha
                rgba = np.asarray(img.convert("RGBA"))
        except Exception as e:
            print(f"Erreur de lecture sur {nom_fichier} : {e}")
            continue

        h, w, _ = rgba.shape
        if h < 450 or w < 950:
            print(f"Image ignorée (trop petite) : {nom_fichier} ({w}x{h})")
            continue

        # Séparation du masque d'opacité et de la végétation
        # Canal 3 = Alpha (opaque si > 127), Canal 0 = R/Luminance (blanc si > 127)
        matrice_valide = rgba[:, :, 3] > 127
        matrice_blanche = (rgba[:, :, 0] > 127) & matrice_valide

        # Découpe des 3 carrés pour les deux masques
        carres_blancs = [
            matrice_blanche[150:450, 50:350],
            matrice_blanche[150:450, 350:650],
            matrice_blanche[150:450, 650:950]
        ]
        carres_valides = [
            matrice_valide[150:450, 50:350],
            matrice_valide[150:450, 350:650],
            matrice_valide[150:450, 650:950]
        ]

        stats = [analyser_carre(cb, cv) for cb, cv in zip(carres_blancs, carres_valides)]

        ligne = [nom_fichier, date_formatee, bac]

        # Données de chaque carré
        for s in stats:
            for cle in CLES_METRIQUES:
                ligne.append(round(s[cle], 4))

        # Moyenne et écart-type inter-carrés
        for cle in CLES_METRIQUES:
            valeurs = [s[cle] for s in stats]
            ligne.append(round(float(np.mean(valeurs)), 4))
            ligne.append(round(float(np.std(valeurs, ddof=1)), 4))

        donnees_csv.append(ligne)

    with open(FICHIER_CSV, mode="w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(en_tetes)
        writer.writerows(donnees_csv)

    print(f"Terminé : {len(donnees_csv)} photo(s) analysée(s) dans {FICHIER_CSV}")

if __name__ == "__main__":
    traiter_images_dossier()