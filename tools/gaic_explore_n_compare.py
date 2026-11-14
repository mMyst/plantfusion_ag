import os
import glob
import re
import math
from datetime import datetime, timedelta
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from matplotlib.patches import Circle
from matplotlib.ticker import FuncFormatter
from openalea.plantgl.all import Scene, TriangleSet


def extract_mesh_data(file_path):
    """Extrait les barycentres XY et les surfaces 3D des triangles d'un fichier."""
    scene = Scene(file_path)
    all_barycenters = []
    all_areas = []

    for shape in scene:
        geom = shape.geometry
        if not isinstance(geom, TriangleSet):
            continue

        pts = np.array([[p.x, p.y, p.z] for p in geom.pointList])
        indices = np.array([[tri[0], tri[1], tri[2]] for tri in geom.indexList])

        if len(indices) == 0:
            continue

        v0, v1, v2 = pts[indices[:, 0]], pts[indices[:, 1]], pts[indices[:, 2]]
        
        barys = (v0 + v1 + v2) / 3.0
        
        cross_products = np.cross(v1 - v0, v2 - v0)
        areas = 0.5 * np.linalg.norm(cross_products, axis=1)

        all_barycenters.append(barys[:, :2])
        all_areas.append(areas)

    if all_barycenters:
        return np.vstack(all_barycenters), np.concatenate(all_areas)
    else:
        return np.array([]), np.array([])


def load_lubbac_mapping(filepath):
    """Crée un dictionnaire associant l'heure t au 'jour continu'."""
    try:
        df_lubbac = pd.read_csv(filepath, sep=None, engine='python')
        base_year = df_lubbac['AN'].min()
        cday = df_lubbac['DOY'] + (df_lubbac['AN'] - base_year) * 365 + (df_lubbac['hour'] / 24.0)
        return dict(zip(df_lubbac['t'], cday))
    except Exception as e:
        print(f"Erreur lors du chargement de LUBBAC: {e}")
        return {}


def load_and_sync_datasets(base_folders_dict, lubbac_mapping):
    """Charge les bgeom et les intervalles foliaires Darwinkel."""
    dataset = {}
    for series_name, base_folder in base_folders_dict.items():
        if not os.path.exists(base_folder):
            print(f"Attention: Le dossier '{base_folder}' n'existe pas.")
            continue
            
        bgeom_folder = os.path.join(base_folder, 'light', 'plantgl')
        axes_file = os.path.join(base_folder, 'wheat', 'postprocessing', 'axes_postprocessing.csv')
        
        # --- 1. Chargement des géométries ---
        files = glob.glob(os.path.join(bgeom_folder, '*.bgeom'))
        parsed_files = []
        for f in files:
            basename = os.path.basename(f)
            match = re.search(r'(\d+)(\d{3})(\d{2})[^\d]*\.bgeom$', basename.lower())
            if match:
                xxx, ddd, hh = int(match.group(1)), int(match.group(2)), int(match.group(3))
                parsed_files.append((xxx, ddd, hh, f, basename))
        
        parsed_files.sort(key=lambda item: (item[0], item[2]))
        
        folder_data = []
        seen_days = set()
        year_offset = 0
        prev_ddd = None
        
        for xxx, ddd, hh, filepath, basename in parsed_files:
            if prev_ddd is not None and ddd < prev_ddd - 100:
                year_offset += 365
                
            continuous_day = ddd + year_offset + (hh / 24.0)
            prev_ddd = ddd
            
            day_key = ddd + year_offset
            if day_key not in seen_days:
                seen_days.add(day_key)
                barys, areas = extract_mesh_data(filepath)
                folder_data.append({
                    'filename': basename,
                    'continuous_day': continuous_day,
                    'barys': barys,
                    'areas': areas
                })
                
        # --- 2. Chargement des intervalles foliaires ---
        series_intervals = []
        if os.path.exists(axes_file):
            try:
                df_axes = pd.read_csv(axes_file)
                if 't' in df_axes.columns and 'nb_leaves' in df_axes.columns:
                    df_sum = df_axes.groupby('t')['nb_leaves'].sum().reset_index().sort_values('t')
                    change_mask = df_sum['nb_leaves'].diff() != 0
                    change_mask.iloc[0] = True
                    df_changes = df_sum[change_mask]
                    
                    t_max = df_sum['t'].max()
                    for i in range(len(df_changes)):
                        t_start = df_changes.iloc[i]['t']
                        val = df_changes.iloc[i]['nb_leaves']
                        t_end = df_changes.iloc[i+1]['t'] if i+1 < len(df_changes) else t_max
                        
                        if t_start in lubbac_mapping and t_end in lubbac_mapping:
                            c_start = lubbac_mapping[t_start]
                            c_end = lubbac_mapping[t_end]
                            if c_end > c_start:
                                series_intervals.append({'start': c_start, 'end': c_end, 'val': val})
            except Exception as e:
                print(f"Erreur de lecture de {axes_file} : {e}")
        else:
            print(f"Attention: Fichier axes_postprocessing introuvable pour {series_name}.")

        dataset[series_name] = {'bgeom': folder_data, 'intervals': series_intervals}
        print(f"[{series_name}] : {len(folder_data)} jours et {len(series_intervals)} intervalles chargés.")
        
    return dataset


def compute_metrics_series(bgeom_data, target_x, target_y, radius):
    """Calcule instantanément les métriques GAIp pour l'axe continu."""
    days, gaips, total_areas, counts = [], [], [], []
    target = np.array([target_x, target_y])
    circle_area = math.pi * (radius ** 2)

    for data in bgeom_data:
        barys = data['barys']
        areas = data['areas']
        
        if len(barys) > 0:
            distances = np.linalg.norm(barys - target, axis=1)
            mask = distances <= radius
            tot_area, count = np.sum(areas[mask]), np.sum(mask)
            gaip = tot_area / circle_area if circle_area > 0 else 0.0
        else:
            tot_area, count, gaip = 0.0, 0, 0.0
            
        days.append(data['continuous_day'])
        gaips.append(gaip)
        total_areas.append(tot_area)
        counts.append(count)
        
    return days, gaips, total_areas, counts


def format_day_to_date(x, pos):
    day_of_year = int(x) % 365
    if day_of_year == 0 and int(x) > 0: 
        day_of_year = 365
    dt = datetime(2023, 1, 1) + timedelta(days=day_of_year - 1)
    return dt.strftime("%d/%m")


# ==========================================
# CONFIGURATION PRINCIPALE
# ==========================================
if __name__ == "__main__":
    
    # 1. Chemin vers LUBBAC
    LUBBAC_PATH = r"C:\Users\agrumel\Code\Python_Ecophy\plantfusion_ag\inputs_fspmwheat\LUBBAC_H_24_25_capteurs.csv"
    
    # 2. Dossiers de BASE des séries
    FOLDERS_BASE = {
        "Wheat": r"C:\Users\agrumel\Code\Python_Ecophy\plantfusion_ag\outputs\wheat_LUBBAC_emergence_local", 
        "WheatLuz": r"C:\Users\agrumel\Code\Python_Ecophy\plantfusion_ag\outputs\WheatLuz_LUBBAC_emergence_local",
     }

    FOLDERS_BASE = {
        "Wheat": r"C:\Users\agrumel\Documents\Data\Simul LUBBAC\Wheat_LUBBAC_emergence_adel", 
        "WheatLuz": r"C:\Users\agrumel\Documents\Data\Simul LUBBAC\WheatLuz_LUBBAC_emergence_adel",
        "WheatLuzREG": r"C:\Users\agrumel\Documents\Data\Simul LUBBAC\WheatLuzREG_LUBBAC_emergence_adel_autoGAI",
     }
    # 3. Périodes temporelles d'intérêt (DOY)
    PERIODS = {
        'T1': (329, 332), # 25/11 au 28/11
        'T2': (337, 342), # 03/12 au 08/12
        'T3': (343, 352), # 09/12 au 18/12
        'T4': (358, 365),  # 24/12 au 31/12
        'T5': (365+6, 365+11)  # 06/01 au 11/01
    }
    
    INIT_TARGET_X, INIT_TARGET_Y, INIT_RADIUS, INIT_GAIC = 0.125, 0.125, 0.10, 0.07

    lubbac_mapping = load_lubbac_mapping(LUBBAC_PATH)
    dataset = load_and_sync_datasets(FOLDERS_BASE, lubbac_mapping)
    
    if not dataset:
        exit("Aucune donnée chargée. Vérifiez les chemins.")

    fig = plt.figure(figsize=(15, 9))
    fig.canvas.manager.set_window_title('Analyse GAIp & Dynamique Foliaire')
    
    gs = fig.add_gridspec(2, 3, bottom=0.30, hspace=0.3, wspace=0.3)
    ax_info = fig.add_subplot(gs[0, 0])
    ax_info.axis('off')
    text_crossover = ax_info.text(0.5, 0.5, '', transform=ax_info.transAxes, 
                                 horizontalalignment='center', verticalalignment='center', fontsize=11,
                                 bbox=dict(boxstyle='round,pad=1', facecolor='#f8f9fa', alpha=1.0, edgecolor='gray'))

    ax_visu = fig.add_subplot(gs[1, 0])
    ax_gai = fig.add_subplot(gs[0, 1:])
    ax_area = fig.add_subplot(gs[1, 1])
    ax_count = fig.add_subplot(gs[1, 2])

    # --- Schéma de la scène ---
    ax_visu.set_xlim(0, 0.30)
    ax_visu.set_ylim(0, 0.30)
    ax_visu.set_xticks(np.arange(0, 0.35, 0.05))
    ax_visu.set_yticks(np.arange(0, 0.35, 0.05))
    ax_visu.grid(True, linestyle='--', color='gray', alpha=0.7)
    ax_visu.set_title("Scène vue de haut (30x30 cm)")
    ax_visu.set_aspect('equal')

    target_point, = ax_visu.plot(INIT_TARGET_X, INIT_TARGET_Y, 'ro', markersize=6)
    radius_circle = Circle((INIT_TARGET_X, INIT_TARGET_Y), INIT_RADIUS, fill=True, color='red', alpha=0.2)
    ax_visu.add_patch(radius_circle)

    START_DAY_CONTINUOUS = 322 # 18 Novembre
    END_DAY_CONTINUOUS = 395   # 30 Janvier

    # --- Configuration des courbes et INTERVALLES FOLIAIRES ---
    lines_gai, lines_area, lines_count = {}, {}, {}
    colors = plt.cm.tab10.colors
    
    for idx, (series_name, data_pack) in enumerate(dataset.items()):
        c = colors[idx % len(colors)]
        lines_gai[series_name], = ax_gai.plot([], [], marker='o', markersize=4, color=c, label=series_name)
        lines_area[series_name], = ax_area.plot([], [], marker='s', markersize=4, color=c)
        lines_count[series_name], = ax_count.plot([], [], marker='^', markersize=4, color=c)

        y_text_base = 0.98 - (idx * 0.08) 
        
        for inter in data_pack['intervals']:
            ax_gai.axvline(x=inter['start'], color=c, linestyle='--', alpha=0.3, linewidth=1)
            ax_gai.axvspan(inter['start'], inter['end'], color=c, alpha=0.05)
            
            vis_start = max(inter['start'], START_DAY_CONTINUOUS)
            vis_end = min(inter['end'], END_DAY_CONTINUOUS)
            center_x = (inter['start'] + inter['end']) / 2 if vis_start >= vis_end else (vis_start + vis_end) / 2
                
            label_text = f"F{int(inter['val'])}"
            ax_gai.text(center_x, y_text_base, label_text,
                        color=c, transform=ax_gai.get_xaxis_transform(),
                        ha='center', va='top', fontsize=8, alpha=0.9, fontweight='bold',
                        bbox=dict(facecolor='#f8f9fa', edgecolor=c, alpha=0.7, boxstyle='round,pad=0.2', linewidth=0.5))

    # --- AJOUT DES PÉRIODES T1 à T4 ---
    for period_name, (start_doy, end_doy) in PERIODS.items():
        # Applique l'ombrage sur les 3 graphiques
        for ax in [ax_gai, ax_area, ax_count]:
            ax.axvspan(start_doy, end_doy, color='#FFD700', alpha=0.15, zorder=0)
        
        # Ajoute l'étiquette T1..T4 en bas du graphique principal
        ax_gai.text((start_doy + end_doy) / 2, 0.05, period_name, 
                    transform=ax_gai.get_xaxis_transform(),
                    ha='center', va='bottom', fontsize=9, fontweight='bold', color='#555555',
                    bbox=dict(facecolor='white', edgecolor='gray', alpha=0.8, boxstyle='round,pad=0.2', linewidth=0.5))

    hline_gaic = ax_gai.axhline(INIT_GAIC, color='red', linestyle='--', alpha=0.8, label="Seuil GAIc")
    date_formatter = FuncFormatter(format_day_to_date)
    
    for ax in [ax_gai, ax_area, ax_count]:
        ax.xaxis.set_major_formatter(date_formatter)
        ax.set_xlim(START_DAY_CONTINUOUS, END_DAY_CONTINUOUS)
        ax.grid(True, alpha=0.3)
        ax.set_xlabel("Date (JJ/MM)")

    ax_gai.set_ylim(0, 0.30)
    ax_area.set_ylim(0, 0.04)
    ax_count.set_ylim(0, 2000)

    ax_gai.set_title("Évolution du GAIp")
    ax_gai.legend()
    ax_area.set_title("Surface des triangles ciblés (m²)")
    ax_count.set_title("Nombre de triangles")

    # --- Création des Sliders ---
    ax_slider_c = fig.add_axes([0.15, 0.20, 0.65, 0.03])
    ax_slider_x = fig.add_axes([0.15, 0.15, 0.65, 0.03])
    ax_slider_y = fig.add_axes([0.15, 0.10, 0.65, 0.03])
    ax_slider_r = fig.add_axes([0.15, 0.05, 0.65, 0.03])

    slider_c = Slider(ax_slider_c, 'Seuil GAIc', 0.0, 0.50, valinit=INIT_GAIC, valstep=0.01)
    slider_x = Slider(ax_slider_x, 'Target X (m)', 0.0, 0.30, valinit=INIT_TARGET_X, valstep=0.005)
    slider_y = Slider(ax_slider_y, 'Target Y (m)', 0.0, 0.30, valinit=INIT_TARGET_Y, valstep=0.005)
    slider_r = Slider(ax_slider_r, 'Radius (m)', 0.01, 0.20, valinit=INIT_RADIUS, valstep=0.005)

    def update(val):
        tx, ty, r = slider_x.val, slider_y.val, slider_r.val
        gaic_threshold = slider_c.val

        target_point.set_data([tx], [ty])
        radius_circle.center = (tx, ty)
        radius_circle.set_radius(r)
        
        hline_gaic.set_ydata([gaic_threshold, gaic_threshold])

        crossover_info = [f"** Dépassement du seuil (GAIp > {gaic_threshold:.2f}) **\n"]

        for series_name, data_pack in dataset.items():
            days, gaips, areas, counts = compute_metrics_series(data_pack['bgeom'], tx, ty, r)
            
            lines_gai[series_name].set_data(days, gaips)
            lines_area[series_name].set_data(days, areas)
            lines_count[series_name].set_data(days, counts)

            gaips_arr = np.array(gaips)
            indices = np.where(gaips_arr > gaic_threshold)[0]
            
            if len(indices) > 0:
                first_date = format_day_to_date(days[indices[0]], None)
                crossover_info.append(f"• {series_name} : {first_date}")
            else:
                crossover_info.append(f"• {series_name} : Non atteint")

        text_crossover.set_text("\n".join(crossover_info))
        fig.canvas.draw_idle()

    slider_c.on_changed(update)
    slider_x.on_changed(update)
    slider_y.on_changed(update)
    slider_r.on_changed(update)

    update(None)
    plt.show()

