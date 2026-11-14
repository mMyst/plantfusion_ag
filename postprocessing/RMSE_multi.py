"""
RMSE_multi.py
-------------
Script d'évaluation et de comparaison multi-modèles basé sur le modèle de référence 'Control'.
Remplace les courbes temporelles par des graphiques de parité (1:1 / Scatter plots) 
avec calcul automatique du R², du RMSE et du nRMSE (%) : nRMSE = (RMSE / y_mean_ctrl) * 100.

Auteur: Adapté à partir de comparison_multi.py
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


# ==============================================================================
# Fonctions de calcul des métriques statistiques (R², RMSE, nRMSE)
# ==============================================================================

def compute_metrics(y_control, y_model):
    """
    Calcule le R², le RMSE et le nRMSE (%) entre le Control (référence/observé) et le modèle simulé.
    Formule nRMSE : nRMSE = (RMSE / y_mean_control) * 100
    
    Paramètres:
        y_control (array-like): Valeurs du modèle Control (référence / valeurs observées).
        y_model (array-like): Valeurs du modèle simulé à évaluer.
        
    Retourne:
        tuple: (r2, rmse, nrmse)
    """
    yc = np.asarray(y_control, dtype=float)
    ym = np.asarray(y_model, dtype=float)
    
    # Masque pour exclure les valeurs manquantes ou infinies
    mask = np.isfinite(yc) & np.isfinite(ym)
    yc = yc[mask]
    ym = ym[mask]
    
    if len(yc) < 2:
        return np.nan, np.nan, np.nan
    
    # Root Mean Square Error (RMSE)
    rmse = np.sqrt(np.mean((ym - yc) ** 2))
    
    # Moyenne de la référence Control (y_bar) et calcul du nRMSE (%)
    mean_c = np.mean(yc)
    if abs(mean_c) > 1e-12:
        nrmse = (rmse / mean_c) * 100.0
    else:
        nrmse = np.nan
    
    # Coefficient de détermination R² (Pearson r au carré)
    var_c = np.var(yc)
    var_m = np.var(ym)
    if var_c > 1e-12 and var_m > 1e-12:
        corr_matrix = np.corrcoef(yc, ym)
        r = corr_matrix[0, 1]
        r2 = r ** 2
    else:
        r2 = np.nan
        
    return r2, rmse, nrmse


# Liste globale pour collecter toutes les métriques et générer le tableau récapitulatif
METRICS_RECORDS = []

def record_metric(variable, sim_name, r2, rmse, nrmse, category=""):
    """Enregistre les métriques calculées pour le rapport de synthèse."""
    METRICS_RECORDS.append({
        'Category': category,
        'Variable': variable,
        'Simulation': sim_name,
        'R2': r2,
        'RMSE': rmse,
        'nRMSE (%)': nrmse
    })


# ==============================================================================
# Fonction générique de tracé de parité (1:1) avec R², RMSE et nRMSE
# ==============================================================================

def plot_parity(ax, df_control, dict_sims, x_col, y_col=None, merge_on='t', 
                label_x='Control', var_title='', unit='', category=''):
    """
    Trace un graphique 1:1 pour une variable donnée en comparant plusieurs modèles au Control.
    Affiche R², RMSE et nRMSE en légende pour chaque modèle.
    """
    if y_col is None:
        y_col = x_col

    all_vals_c = []
    all_vals_m = []

    # Palette de couleurs automatique
    colors = plt.cm.tab10.colors

    for idx, (sim_name, df_sim) in enumerate(dict_sims.items()):
        color = colors[idx % len(colors)]
        
        # Jointure interne sur l'identifiant temporel (ou métamère) pour aligner exactement les points
        if isinstance(merge_on, list):
            df_merged = pd.merge(df_control[merge_on + [x_col]], 
                                 df_sim[merge_on + [y_col]], 
                                 on=merge_on, 
                                 suffixes=('_ctrl', '_sim'))
        else:
            df_merged = pd.merge(df_control[[merge_on, x_col]], 
                                 df_sim[[merge_on, y_col]], 
                                 on=merge_on, 
                                 suffixes=('_ctrl', '_sim'))

        col_c = f"{x_col}_ctrl" if x_col == y_col else x_col
        col_m = f"{y_col}_sim" if x_col == y_col else y_col
        
        # Nettoyage des données
        sub = df_merged[[col_c, col_m]].dropna()
        if len(sub) == 0:
            continue
            
        yc = sub[col_c].values
        ym = sub[col_m].values
        
        r2, rmse, nrmse = compute_metrics(yc, ym)
        record_metric(var_title, sim_name, r2, rmse, nrmse, category=category)
        
        # Format d'affichage des métriques
        r2_str = f"{r2:.2f}" if np.isfinite(r2) else "NA"
        if np.isfinite(rmse):
            rmse_str = f"{rmse:.2e}" if (abs(rmse) < 0.01 or abs(rmse) > 1000) else f"{rmse:.2f}"
        else:
            rmse_str = "NA"
        nrmse_str = f"{nrmse:.1f}%" if np.isfinite(nrmse) else "NA"
            
        legend_label = f"{sim_name} (R²={r2_str}, RMSE={rmse_str}, nRMSE={nrmse_str})"
        
        # Tracé des points
        ax.scatter(yc, ym, label=legend_label, color=color, alpha=0.65, s=22, edgecolors='none')
        
        all_vals_c.extend(yc)
        all_vals_m.extend(ym)

    # Tracé de la ligne 1:1 (diagonale de référence parfaite)
    if len(all_vals_c) > 0 and len(all_vals_m) > 0:
        min_val = min(np.min(all_vals_c), np.min(all_vals_m))
        max_val = max(np.max(all_vals_c), np.max(all_vals_m))
        span = max_val - min_val
        if span == 0:
            span = 1.0
        min_plot = min_val - 0.05 * span
        max_plot = max_val + 0.05 * span
        
        ax.plot([min_plot, max_plot], [min_plot, max_plot], 'k--', lw=1.2, label='1:1', zorder=1)
        ax.set_xlim(min_plot, max_plot)
        ax.set_ylim(min_plot, max_plot)

    # Formatage des axes et labels
    unit_str = f" ({unit})" if unit else ""
    ax.set_xlabel(f"{label_x}{unit_str}", fontsize=9)
    ax.set_ylabel(f"Simulations{unit_str}", fontsize=9)
    ax.set_title(var_title, fontsize=10, fontweight='bold')
    ax.grid(True, linestyle=':', alpha=0.6)
    ax.legend(fontsize=6.5, loc='best', framealpha=0.85)


# ==============================================================================
# Fonctions de modules thématiques (Phloem, Roots, Dry Mass, etc.)
# ==============================================================================

def phloem(dict_current_organs, df_control_organs, pdf):
    fig, axs = plt.subplots(2, 2, figsize=(11, 9))

    df_c_phloem = df_control_organs[df_control_organs.organ == 'phloem']
    dict_sim_phloem = {name: df[df.organ == 'phloem'] for name, df in dict_current_organs.items()}

    # 1. Concentration sucrose
    plot_parity(axs[0, 0], df_c_phloem, dict_sim_phloem, 
                x_col='Conc_Sucrose', var_title='Phloem Conc_Sucrose', 
                unit='µmol g-1', category='Phloem')

    # 2. Amount of sucrose
    plot_parity(axs[0, 1], df_c_phloem, dict_sim_phloem, 
                x_col='sucrose', var_title='Phloem Amount Sucrose', 
                unit='µmol C', category='Phloem')

    # 3. Concentration amino acids
    plot_parity(axs[1, 0], df_c_phloem, dict_sim_phloem, 
                x_col='Conc_Amino_Acids', var_title='Phloem Conc_Amino_Acids', 
                unit='µmol g-1', category='Phloem')

    # 4. Amino acids amount
    plot_parity(axs[1, 1], df_c_phloem, dict_sim_phloem, 
                x_col='amino_acids', var_title='Phloem Amount Amino Acids', 
                unit='µmol N', category='Phloem')

    fig.suptitle("Phloème : R², RMSE et nRMSE (Modèles vs Control)", fontsize=13, fontweight='bold')
    plt.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def photosynthesis(dict_current_axes, df_control_axes, pdf):
    fig, axs = plt.subplots(1, 2, figsize=(12, 5))

    # Agrégation journalière pour Control
    df_c = df_control_axes.copy()
    df_c['day'] = df_c['t'] // 24 + 1
    df_c_daily = df_c.groupby('day')[['Total_Photosynthesis', 'Tillers_Photosynthesis']].sum().reset_index()

    # Agrégation journalière pour les simulations
    dict_sim_daily = {}
    for name, df in dict_current_axes.items():
        df_sim = df.copy()
        df_sim['day'] = df_sim['t'] // 24 + 1
        dict_sim_daily[name] = df_sim.groupby('day')[['Total_Photosynthesis', 'Tillers_Photosynthesis']].sum().reset_index()

    # 1. Total Photosynthesis
    plot_parity(axs[0], df_c_daily, dict_sim_daily, 
                x_col='Total_Photosynthesis', merge_on='day', label_x='Control (day sum)',
                var_title='Total Photosynthesis', unit='µmol C', category='Photosynthesis')

    # 2. Tillers Photosynthesis
    plot_parity(axs[1], df_c_daily, dict_sim_daily, 
                x_col='Tillers_Photosynthesis', merge_on='day', label_x='Control (day sum)',
                var_title='Tillers Photosynthesis', unit='µmol C', category='Photosynthesis')

    fig.suptitle("Photosynthèse : R², RMSE et nRMSE (Modèles vs Control)", fontsize=13, fontweight='bold')
    plt.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def roots(dict_current_organs, df_control_organs, pdf):
    fig, axs = plt.subplots(2, 2, figsize=(11, 9))

    df_c_roots = df_control_organs[df_control_organs.organ == 'roots']
    dict_sim_roots = {name: df[df.organ == 'roots'] for name, df in dict_current_organs.items()}

    # 1. Concentration sucrose
    plot_parity(axs[0, 0], df_c_roots, dict_sim_roots, 
                x_col='Conc_Sucrose', var_title='Roots Conc_Sucrose', 
                unit='µmol g-1', category='Roots')

    # 2. Amount of sucrose
    plot_parity(axs[0, 1], df_c_roots, dict_sim_roots, 
                x_col='sucrose', var_title='Roots Amount Sucrose', 
                unit='µmol C', category='Roots')

    # 3. Concentration nitrates
    plot_parity(axs[1, 0], df_c_roots, dict_sim_roots, 
                x_col='Conc_Nitrates', var_title='Roots Conc_Nitrates', 
                unit='µmol g-1', category='Roots')

    # 4. Cytokinins
    plot_parity(axs[1, 1], df_c_roots, dict_sim_roots, 
                x_col='Conc_cytokinins', var_title='Roots Conc_cytokinins', 
                unit='AU g-1', category='Roots')

    fig.suptitle("Racines : R², RMSE et nRMSE (Modèles vs Control)", fontsize=13, fontweight='bold')
    plt.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def dry_mass(dict_current_axes, df_control_axes, dict_current_organs, df_control_organs, pdf):
    fig, axs = plt.subplots(2, 2, figsize=(11, 9))

    # 1. Dry mass shoot
    plot_parity(axs[0, 0], df_control_axes, dict_current_axes, 
                x_col='sum_dry_mass_shoot', var_title='Dry Mass Shoot', 
                unit='g', category='Biomass')

    # 2. Dry mass roots
    plot_parity(axs[0, 1], df_control_axes, dict_current_axes, 
                x_col='sum_dry_mass_roots', var_title='Dry Mass Roots', 
                unit='g', category='Biomass')

    # 3. mstruct shoot
    plot_parity(axs[1, 0], df_control_axes, dict_current_axes, 
                x_col='mstruct_shoot', var_title='mstruct Shoot', 
                unit='g', category='Biomass')

    # 4. mstruct roots
    df_c_roots = df_control_organs[df_control_organs['organ'] == 'roots']
    dict_sim_roots = {name: df[df['organ'] == 'roots'] for name, df in dict_current_organs.items()}
    plot_parity(axs[1, 1], df_c_roots, dict_sim_roots, 
                x_col='mstruct', var_title='mstruct Roots', 
                unit='g', category='Biomass')

    fig.suptitle("Biomasses sèches : R², RMSE et nRMSE (Modèles vs Control)", fontsize=13, fontweight='bold')
    plt.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)

    # 5. Ratio shoot:root
    fig, ax = plt.subplots(figsize=(7, 6))
    plot_parity(ax, df_control_axes, dict_current_axes, 
                x_col='shoot_roots_ratio', var_title='Shoot : Root Ratio', 
                unit='', category='Biomass')
    fig.suptitle("Ratio Shoot / Root : R², RMSE et nRMSE", fontsize=12, fontweight='bold')
    plt.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def N_mass(dict_current_axes, df_control_axes, dict_current_organs, df_control_organs, pdf):
    fig, axs = plt.subplots(2, 2, figsize=(11, 9))

    # 1. % N axis
    plot_parity(axs[0, 0], df_control_axes, dict_current_axes, 
                x_col='N_content', var_title='N Content Axis', 
                unit='% DM', category='Nitrogen')

    # 2. N shoot
    plot_parity(axs[0, 1], df_control_axes, dict_current_axes, 
                x_col='N_content_shoot', var_title='N Content Shoot', 
                unit='% DM', category='Nitrogen')

    # 3. N axis (g)
    plot_parity(axs[1, 0], df_control_axes, dict_current_axes, 
                x_col='sum_N_g', var_title='Total N Axis', 
                unit='g', category='Nitrogen')

    # 4. N uptake roots
    df_c_roots = df_control_organs[df_control_organs['organ'] == 'roots']
    dict_sim_roots = {name: df[df['organ'] == 'roots'] for name, df in dict_current_organs.items()}
    plot_parity(axs[1, 1], df_c_roots, dict_sim_roots, 
                x_col='Uptake_Nitrates', var_title='Nitrate Uptake Roots', 
                unit='µmol', category='Nitrogen')

    fig.suptitle("Azote (N) : R², RMSE et nRMSE (Modèles vs Control)", fontsize=13, fontweight='bold')
    plt.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def surface(dict_current_elements, df_control_elements, pdf):
    fig, axs = plt.subplots(2, 2, figsize=(11, 9))

    # Total green area
    df_c_tot = df_control_elements.groupby('t')['green_area'].sum().reset_index()
    dict_sim_tot = {name: df.groupby('t')['green_area'].sum().reset_index() for name, df in dict_current_elements.items()}
    plot_parity(axs[0, 0], df_c_tot, dict_sim_tot, x_col='green_area', 
                var_title='Total Green Area', unit='m²', category='Surface')

    # Blade green area
    df_c_b = df_control_elements[df_control_elements.organ == 'blade'].groupby('t')['green_area'].sum().reset_index()
    dict_sim_b = {name: df[df.organ == 'blade'].groupby('t')['green_area'].sum().reset_index() for name, df in dict_current_elements.items()}
    plot_parity(axs[0, 1], df_c_b, dict_sim_b, x_col='green_area', 
                var_title='Blade Green Area', unit='m²', category='Surface')

    # Sheath green area
    df_c_s = df_control_elements[df_control_elements.organ == 'sheath'].groupby('t')['green_area'].sum().reset_index()
    dict_sim_s = {name: df[df.organ == 'sheath'].groupby('t')['green_area'].sum().reset_index() for name, df in dict_current_elements.items()}
    plot_parity(axs[1, 0], df_c_s, dict_sim_s, x_col='green_area', 
                var_title='Sheath Green Area', unit='m²', category='Surface')

    # Internode green area
    df_c_i = df_control_elements[df_control_elements.organ == 'internode'].groupby('t')['green_area'].sum().reset_index()
    dict_sim_i = {name: df[df.organ == 'internode'].groupby('t')['green_area'].sum().reset_index() for name, df in dict_current_elements.items()}
    plot_parity(axs[1, 1], df_c_i, dict_sim_i, x_col='green_area', 
                var_title='Internode Green Area', unit='m²', category='Surface')

    fig.suptitle("Surfaces vertes foliaires : R², RMSE et nRMSE (Modèles vs Control)", fontsize=13, fontweight='bold')
    plt.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def height(dict_current_elements_outputs, df_control_elements_outputs, pdf):
    fig, ax = plt.subplots(figsize=(7, 6))

    df_c_h = df_control_elements_outputs.groupby('t')['height'].max().reset_index()
    dict_sim_h = {name: df.groupby('t')['height'].max().reset_index() for name, df in dict_current_elements_outputs.items()}

    plot_parity(ax, df_c_h, dict_sim_h, x_col='height', 
                var_title='Plant Height', unit='m', category='Morphology')

    fig.suptitle("Hauteur de la plante : R², RMSE et nRMSE", fontsize=12, fontweight='bold')
    plt.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def plastochrone(dict_current_axes_outputs, df_control_axes_outputs, pdf):
    fig, ax = plt.subplots(figsize=(7, 6))

    plot_parity(ax, df_control_axes_outputs, dict_current_axes_outputs, 
                x_col='nb_leaves', var_title='Number of Leaves on MS', 
                unit='leaves', category='Development')

    fig.suptitle("Plastochrone (Nombre de feuilles) : R², RMSE et nRMSE", fontsize=12, fontweight='bold')
    plt.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def leaf_emergence(dict_current_hz, df_control_hz, pdf, delta_t=0):
    """
    Compare les temps d'émergence des feuilles entre le Control et chaque modèle.
    """
    t_emergence_ctrl = {}
    for phyto_id in df_control_hz['metamer'].unique():
        if phyto_id == 0 or not df_control_hz[(df_control_hz.metamer == phyto_id)]['leaf_is_emerged'].any():
            continue
        sub = df_control_hz[(df_control_hz.metamer == phyto_id) & (df_control_hz.leaf_is_emerged == True)]
        t_emergence_ctrl[phyto_id] = sub['t'].iloc[0] - delta_t

    df_ctrl_em = pd.DataFrame(list(t_emergence_ctrl.items()), columns=['metamer', 't_emergence'])

    dict_sim_em = {}
    for name, df in dict_current_hz.items():
        t_em = {}
        for phyto_id in df['metamer'].unique():
            if phyto_id == 0 or not df[(df.metamer == phyto_id)]['leaf_is_emerged'].any():
                continue
            sub = df[(df.metamer == phyto_id) & (df.leaf_is_emerged == True)]
            t_em[phyto_id] = sub['t'].iloc[0] - delta_t
        dict_sim_em[name] = pd.DataFrame(list(t_em.items()), columns=['metamer', 't_emergence'])

    fig, ax = plt.subplots(figsize=(7, 6))
    plot_parity(ax, df_ctrl_em, dict_sim_em, x_col='t_emergence', merge_on='metamer',
                label_x='Control Emergence', var_title='Leaf Emergence Time', 
                unit='h', category='Development')

    fig.suptitle("Temps d'émergence foliaire : R², RMSE et nRMSE", fontsize=12, fontweight='bold')
    plt.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def leaf_length_mstruct_area(dict_current_hz, df_control_hz, dict_current_elements, df_control_elements, pdf):
    """
    Trace pour chaque phytomère/métamère les graphiques de parité 1:1 avec R², RMSE et nRMSE.
    """
    for phyto_id in df_control_hz['metamer'].unique():
        fig, axs = plt.subplots(3, 3, figsize=(13, 11))

        # 1. Length leaf
        if phyto_id in (1, 2):
            df_c_l = df_control_elements[df_control_elements.metamer == phyto_id].groupby('t')['length'].sum().reset_index()
            dict_sim_l = {name: df[df.metamer == phyto_id].groupby('t')['length'].sum().reset_index() 
                          for name, df in dict_current_elements.items()}
            plot_parity(axs[0, 0], df_c_l, dict_sim_l, x_col='length', 
                        var_title=f'Leaf_{phyto_id} Length', unit='m', category=f'Phytomer_{phyto_id}')
        else:
            df_c_l = df_control_hz[df_control_hz.metamer == phyto_id].groupby('t')['leaf_L'].sum().reset_index()
            dict_sim_l = {name: df[df.metamer == phyto_id].groupby('t')['leaf_L'].sum().reset_index() 
                          for name, df in dict_current_hz.items()}
            plot_parity(axs[0, 0], df_c_l, dict_sim_l, x_col='leaf_L', 
                        var_title=f'Leaf_{phyto_id} Length', unit='m', category=f'Phytomer_{phyto_id}')

        # 2. Length sheath
        df_c_sh_l = df_control_elements[(df_control_elements.metamer == phyto_id) & (df_control_elements.organ == 'sheath')].groupby('t')['length'].sum().reset_index()
        dict_sim_sh_l = {name: df[(df.metamer == phyto_id) & (df.organ == 'sheath')].groupby('t')['length'].sum().reset_index() 
                         for name, df in dict_current_elements.items()}
        plot_parity(axs[0, 1], df_c_sh_l, dict_sim_sh_l, x_col='length', 
                    var_title=f'Sheath_{phyto_id} Length', unit='m', category=f'Phytomer_{phyto_id}')

        # 3. Internode length
        df_c_in_l = df_control_hz[df_control_hz.metamer == phyto_id].groupby('t')['internode_L'].sum().reset_index()
        dict_sim_in_l = {name: df[df.metamer == phyto_id].groupby('t')['internode_L'].sum().reset_index() 
                         for name, df in dict_current_hz.items()}
        plot_parity(axs[0, 2], df_c_in_l, dict_sim_in_l, x_col='internode_L', 
                    var_title=f'Internode_{phyto_id} Length', unit='m', category=f'Phytomer_{phyto_id}')

        # 4. Mstruct hz
        if phyto_id in (1, 2):
            df_c_ms = df_control_elements[(df_control_elements.metamer == phyto_id) & (df_control_elements.element == 'LeafElement1')].groupby('t')['mstruct'].sum().reset_index()
            dict_sim_ms = {name: df[(df.metamer == phyto_id) & (df.element == 'LeafElement1')].groupby('t')['mstruct'].sum().reset_index() 
                           for name, df in dict_current_elements.items()}
        else:
            df_c_ms = df_control_hz[df_control_hz.metamer == phyto_id].groupby('t')['mstruct'].sum().reset_index()
            dict_sim_ms = {name: df[df.metamer == phyto_id].groupby('t')['mstruct'].sum().reset_index() 
                           for name, df in dict_current_hz.items()}
        plot_parity(axs[1, 0], df_c_ms, dict_sim_ms, x_col='mstruct', 
                    var_title=f'Mstruct Hz_{phyto_id}', unit='g', category=f'Phytomer_{phyto_id}')

        # 5. mstruct lamina
        if phyto_id != 0:
            df_c_lam = df_control_elements[(df_control_elements.metamer == phyto_id) & (df_control_elements.organ == 'blade')].groupby('t')['mstruct'].sum().reset_index()
            dict_sim_lam = {name: df[(df.metamer == phyto_id) & (df.organ == 'blade')].groupby('t')['mstruct'].sum().reset_index() 
                            for name, df in dict_current_elements.items()}
            plot_parity(axs[1, 1], df_c_lam, dict_sim_lam, x_col='mstruct', 
                        var_title=f'Lamina Mstruct_{phyto_id}', unit='g', category=f'Phytomer_{phyto_id}')
        else:
            axs[1, 1].axis('off')

        # 6. N content lamina
        if phyto_id != 0:
            df_c_n_lam = df_control_elements[(df_control_elements.metamer == phyto_id) & (df_control_elements.element == 'LeafElement1')].groupby('t')['N_tot'].sum().reset_index()
            dict_sim_n_lam = {name: df[(df.metamer == phyto_id) & (df.element == 'LeafElement1')].groupby('t')['N_tot'].sum().reset_index() 
                              for name, df in dict_current_elements.items()}
            plot_parity(axs[1, 2], df_c_n_lam, dict_sim_n_lam, x_col='N_tot', 
                        var_title=f'Lamina N_{phyto_id}', unit='g', category=f'Phytomer_{phyto_id}')
        else:
            axs[1, 2].axis('off')

        # 7. mstruct sheath
        df_c_sh_m = df_control_elements[(df_control_elements.metamer == phyto_id) & (df_control_elements.organ == 'sheath')].groupby('t')['mstruct'].sum().reset_index()
        dict_sim_sh_m = {name: df[(df.metamer == phyto_id) & (df.organ == 'sheath')].groupby('t')['mstruct'].sum().reset_index() 
                         for name, df in dict_current_elements.items()}
        plot_parity(axs[2, 0], df_c_sh_m, dict_sim_sh_m, x_col='mstruct', 
                    var_title=f'Sheath Mstruct_{phyto_id}', unit='g', category=f'Phytomer_{phyto_id}')

        # 8. Green area lamina
        df_c_lam_ga = df_control_elements[(df_control_elements.metamer == phyto_id) & (df_control_elements.element == 'LeafElement1')].groupby('t')['green_area'].sum().reset_index()
        dict_sim_lam_ga = {name: df[(df.metamer == phyto_id) & (df.element == 'LeafElement1')].groupby('t')['green_area'].sum().reset_index() 
                           for name, df in dict_current_elements.items()}
        plot_parity(axs[2, 1], df_c_lam_ga, dict_sim_lam_ga, x_col='green_area', 
                    var_title=f'Lamina GA_{phyto_id}', unit='m²', category=f'Phytomer_{phyto_id}')

        # 9. N content sheath
        if phyto_id != 0:
            df_c_n_sh = df_control_elements[(df_control_elements.metamer == phyto_id) & (df_control_elements.organ == 'sheath')].groupby('t')['N_tot'].sum().reset_index()
            dict_sim_n_sh = {name: df[(df.metamer == phyto_id) & (df.organ == 'sheath')].groupby('t')['N_tot'].sum().reset_index() 
                             for name, df in dict_current_elements.items()}
            plot_parity(axs[2, 2], df_c_n_sh, dict_sim_n_sh, x_col='N_tot', 
                        var_title=f'Sheath N_{phyto_id}', unit='g', category=f'Phytomer_{phyto_id}')
        else:
            axs[2, 2].axis('off')

        fig.suptitle(f"Phytomère {phyto_id} : R², RMSE et nRMSE (Modèles vs Control)", fontsize=13, fontweight='bold')
        plt.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)


def summary_metrics_plots(pdf, output_csv="RMSE_R2_summary.csv"):
    """
    Génère un tableau CSV récapitulatif ainsi que des pages PDF synthétisant
    les scores R², RMSE et nRMSE (%) sous forme de diagrammes en barres groupées.
    """
    if not METRICS_RECORDS:
        return

    df_metrics = pd.DataFrame(METRICS_RECORDS)
    # df_metrics.to_csv(output_csv, index=False)
    # print(f"Tableau récapitulatif des métriques sauvegardé sous : {output_csv}")

    # Synthèse graphique des R², RMSE et nRMSE globaux par catégorie
    categories = df_metrics['Category'].unique()
    for cat in categories:
        if not cat or cat.startswith('Phytomer'):
            continue  # Les phytomères sont très nombreux, on résume les catégories principales
        sub_df = df_metrics[df_metrics['Category'] == cat]
        if sub_df.empty:
            continue

        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 5.5))

        # Pivot pour affichage en barres groupées
        pivot_r2 = sub_df.pivot(index='Variable', columns='Simulation', values='R2')
        pivot_rmse = sub_df.pivot(index='Variable', columns='Simulation', values='RMSE')
        pivot_nrmse = sub_df.pivot(index='Variable', columns='Simulation', values='nRMSE (%)')

        if not pivot_r2.empty:
            pivot_r2.plot(kind='bar', ax=ax1, colormap='tab10', edgecolor='black', alpha=0.85)
            ax1.set_title(f"R² (Corrélation) - {cat}", fontweight='bold')
            ax1.set_ylabel("R²")
            ax1.set_ylim(0, 1.05)
            ax1.grid(True, linestyle=':', alpha=0.6)
            ax1.tick_params(axis='x', rotation=30)
            ax1.legend(fontsize=8)

        if not pivot_rmse.empty:
            pivot_rmse.plot(kind='bar', ax=ax2, colormap='tab10', edgecolor='black', alpha=0.85)
            ax2.set_title(f"RMSE (Erreur brute) - {cat}", fontweight='bold')
            ax2.set_ylabel("RMSE")
            ax2.grid(True, linestyle=':', alpha=0.6)
            ax2.tick_params(axis='x', rotation=30)
            ax2.legend(fontsize=8)

        if not pivot_nrmse.empty:
            pivot_nrmse.plot(kind='bar', ax=ax3, colormap='tab10', edgecolor='black', alpha=0.85)
            ax3.set_title(f"nRMSE (%) (Erreur normalisée) - {cat}", fontweight='bold')
            ax3.set_ylabel("nRMSE (%)")
            ax3.grid(True, linestyle=':', alpha=0.6)
            ax3.tick_params(axis='x', rotation=30)
            ax3.legend(fontsize=8)

        fig.suptitle(f"Synthèse globale des performances : {cat}", fontsize=13, fontweight='bold')
        plt.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)


# ==============================================================================
# Programme principal
# ==============================================================================

if __name__ == '__main__':
    # 1. Dictionnaire des simulations à comparer au Control
    simulations_paths = {
        # 'Default': r'C:\Users\agrumel\Code\Python_Ecophy\plantfusion_ag\outputs\cnwheat_default\wheat',
     
        'Soil3DS_Bound_0.2m': r'C:\Users\agrumel\Code\Python_Ecophy\plantfusion_ag\outputs\cnwheat_soil3ds\bound\0.2m\wheat',
        'Soil3DS_Homogeneous_1.5m': r'C:\Users\agrumel\Code\Python_Ecophy\plantfusion_ag\outputs\cnwheat_soil3ds_debug\homogeneous\1.5m\wheat',
        'Soil3DS_Profile_0.2m': r'C:\Users\agrumel\Code\Python_Ecophy\plantfusion_ag\outputs\cnwheat_soil3ds_debug\profile\0.2m\wheat',
        'Soil3DS_Profile_0.2m_gradient': r'C:\Users\agrumel\Documents\Data\Sorties CNWheat\soil3ds_302010_N/wheat',


        # 'WL_LUBBAC': r'C:\Users\agrumel\Documents\Data\Simul LUBBAC\WheatLuz_LUBBAC_emergence_adel\wheat',
        # 'W_LUBBAC': r'C:\Users\agrumel\Documents\Data\Simul LUBBAC\Wheat_LUBBAC_emergence_adel\wheat',
        
        
    }

    # Chemin vers la simulation de référence 'Control'
    dirpath_control = r'C:\Users\agrumel\Code\Python_Ecophy\WheatFspm\fspm-wheat\example\Vegetative_stages'

    delta_t_simuls = 0
    delta_t_emergence = 0 #1872

    OUTPUTS_CONTROL = os.path.join(dirpath_control, 'brut')
    POSTPROCESSING_CONTROL = os.path.join(dirpath_control, 'postprocessing')

    # Chargement des données de référence (Control)
    print("Chargement des données de référence 'Control'...")
    df_control_axes = pd.read_csv(os.path.join(POSTPROCESSING_CONTROL, 'axes_postprocessing.csv'))
    df_control_axes['t'] = df_control_axes['t'] + delta_t_simuls
    
    df_control_axes_outputs = pd.read_csv(os.path.join(OUTPUTS_CONTROL, 'axes_outputs.csv'))
    df_control_axes_outputs['t'] = df_control_axes_outputs['t'] + delta_t_simuls

    df_control_organs = pd.read_csv(os.path.join(POSTPROCESSING_CONTROL, 'organs_postprocessing.csv'))
    df_control_organs['t'] = df_control_organs['t'] + delta_t_simuls

    df_control_elements = pd.read_csv(os.path.join(POSTPROCESSING_CONTROL, 'elements_postprocessing.csv'))
    df_control_elements['t'] = df_control_elements['t'] + delta_t_simuls

    df_control_elements_outputs = pd.read_csv(os.path.join(OUTPUTS_CONTROL, 'elements_outputs.csv'))
    df_control_elements_outputs['t'] = df_control_elements_outputs['t'] + delta_t_simuls

    df_control_hz = pd.read_csv(os.path.join(POSTPROCESSING_CONTROL, 'hiddenzones_postprocessing.csv'))
    df_control_hz['t'] = df_control_hz['t'] + delta_t_simuls

    # Chargement des simulations à comparer
    print("Chargement des simulations à comparer...")
    dict_current_axes = {}
    dict_current_organs = {}
    dict_current_elements = {}
    dict_current_hz = {}
    dict_current_axes_outputs = {}
    dict_current_elements_outputs = {}

    for sim_name, path in simulations_paths.items():
        post_dir = os.path.join(path, 'postprocessing')
        raw_dir = os.path.join(path, 'brut')

        dict_current_axes[sim_name] = pd.read_csv(os.path.join(post_dir, 'axes_postprocessing.csv'))
        dict_current_organs[sim_name] = pd.read_csv(os.path.join(post_dir, 'organs_postprocessing.csv'))
        dict_current_elements[sim_name] = pd.read_csv(os.path.join(post_dir, 'elements_postprocessing.csv'))
        dict_current_hz[sim_name] = pd.read_csv(os.path.join(post_dir, 'hiddenzones_postprocessing.csv'))

        dict_current_axes_outputs[sim_name] = pd.read_csv(os.path.join(raw_dir, 'axes_outputs.csv'))
        dict_current_elements_outputs[sim_name] = pd.read_csv(os.path.join(raw_dir, 'elements_outputs.csv'))

    # Sous-ensembles pour l'axe principal (MS - Main Stem)
    dict_current_hz_MS = {k: df[df['axis'] == "MS"] for k, df in dict_current_hz.items()}
    df_control_hz_MS = df_control_hz[df_control_hz['axis'] == 'MS']
    
    dict_current_elements_MS = {k: df[df['axis'] == "MS"] for k, df in dict_current_elements.items()}
    df_control_elements_MS = df_control_elements[df_control_elements['axis'] == 'MS']
    
    dict_current_axes_outputs_MS = {k: df[df['axis'] == "MS"] for k, df in dict_current_axes_outputs.items()}
    df_control_axes_outputs_MS = df_control_axes_outputs[df_control_axes_outputs['axis'] == 'MS']

    # 3. Tracé des graphiques de parité R² + RMSE + nRMSE dans un PDF
    pdf_filename = 'RMSE_multi_' + '_'.join(simulations_paths.keys()) + '.pdf'
    print(f"Génération des graphiques R² + RMSE + nRMSE dans : {pdf_filename} ...")
    
    with PdfPages(pdf_filename) as pdf:
        phloem(dict_current_organs, df_control_organs, pdf)
        photosynthesis(dict_current_axes, df_control_axes, pdf)
        roots(dict_current_organs, df_control_organs, pdf)
        dry_mass(dict_current_axes, df_control_axes, dict_current_organs, df_control_organs, pdf)
        N_mass(dict_current_axes, df_control_axes, dict_current_organs, df_control_organs, pdf)
        surface(dict_current_elements, df_control_elements, pdf)
        height(dict_current_elements_outputs, df_control_elements_outputs, pdf)
        leaf_emergence(dict_current_hz_MS, df_control_hz_MS, pdf, delta_t=delta_t_emergence)
        plastochrone(dict_current_axes_outputs_MS, df_control_axes_outputs_MS, pdf)
        leaf_length_mstruct_area(dict_current_hz_MS, df_control_hz_MS, dict_current_elements_MS, df_control_elements_MS, pdf)
        summary_metrics_plots(pdf, output_csv='RMSE_R2_summary.csv')

    print("Terminé avec succès !")