import os
import re
import numpy as np
import pandas as pd
import statsmodels.api as sm
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

# Import optionnel des paramètres elongwheat pour afficher les étoiles "Model parameters"
try:
    from elongwheat import parameters as elongwheat_parameters
    HAS_ELONGWHEAT = True
except ImportError:
    HAS_ELONGWHEAT = False
    print("[Info] Module 'elongwheat' non trouvé : les paramètres théoriques RERmax ne seront pas tracés.")

# ==============================================================================
# 1. CONFIGURATION
# ==============================================================================

# Dictionnaire des simulations à comparer : {'Nom dans la légende': 'Chemin/vers/dossier'}
# Note : Le script cherchera automatiquement dans le dossier ou dans le sous-dossier 'brut/'
SIMULATION_FOLDERS = {
    'Reference': r'C:\Users\agrumel\Code\Python_Ecophy\WheatFspm\fspm-wheat\example\Vegetative_stages',
#    'PF Default': r'C:\Users\agrumel\Code\Python_Ecophy\plantfusion_ag\outputs\cnwheat_default\wheat',
    'Homogeneous': r'C:\Users\agrumel\Code\Python_Ecophy\plantfusion_ag\outputs\cnwheat_soil3ds_debug\homogeneous\1.5m\wheat',
    # 'Bound': r'C:\Users\agrumel\Code\Python_Ecophy\plantfusion_ag\outputs\cnwheat_soil3ds\bound\0.2m\wheat',
    'Profile': r'C:\Users\agrumel\Code\Python_Ecophy\plantfusion_ag\outputs\cnwheat_soil3ds_debug\profile\0.2m\wheat',
    # 'Profile_Ngradient': r'C:\Users\agrumel\Documents\Data\Sorties CNWheat\soil3ds_302010_N\wheat',
    
}

# Chemin vers le fichier d'observations Ljutovac 2002
LJUTOVAC_FILEPATH = r'inputs_fspmwheat\Ljutovac2002.csv'

# Dossier racine des comparaisons + création automatique du sous-dossier à partir des clés
BASE_OUTPUT_DIR = r'C:\Users\agrumel\Documents\Data\Sorties CNWheat\Reference_comparisons'

# Nettoyage des clés pour créer un nom de dossier valide sous Windows/Linux
clean_keys = [re.sub(r'[\\/*?:"<>|\s()]+', '_', k).strip('_') for k in SIMULATION_FOLDERS.keys()]
SUBFOLDER_NAME = "_vs_".join(clean_keys)
OUTPUT_GRAPHS_DIR = os.path.join(BASE_OUTPUT_DIR, SUBFOLDER_NAME)

# Création du sous-dossier de sauvegarde s'il n'existe pas déjà
if not os.path.exists(OUTPUT_GRAPHS_DIR):
    os.makedirs(OUTPUT_GRAPHS_DIR, exist_ok=True)
    print(f"Sous-dossier de sauvegarde créé : '{OUTPUT_GRAPHS_DIR}'")
else:
    print(f"Sous-dossier de sauvegarde déjà existant : '{OUTPUT_GRAPHS_DIR}'")

# Noms des fichiers de sortie bruts nécessaires
HIDDENZONES_FILENAME = "hiddenzones_outputs.csv"
AXES_FILENAME = "axes_outputs.csv"

# Palette de couleurs et marqueurs pour différencier les simulations
COLORS = ['c', 'royalblue', 'forestgreen', 'crimson', 'darkviolet', 'teal', 'magenta', 'brown']
MARKERS = ['o', 's', '^', 'D', 'v', 'p', '*', 'X']


# ==============================================================================
# 2. FONCTIONS DE CALCUL ET DE TRAITEMENT
# ==============================================================================

def resolve_filepath(folder_path, filename):
    """Cherche un fichier dans folder_path/brut/ ou directement dans folder_path/."""
    filepath_brut = os.path.join(folder_path, "brut", filename)
    filepath_direct = os.path.join(folder_path, filename)
    if os.path.isfile(filepath_brut):
        return filepath_brut
    elif os.path.isfile(filepath_direct):
        return filepath_direct
    return None


def compute_metrics(obs, sim):
    """
    Calcule RMSE, nRMSE (par rapport à la moyenne observée) et EF (Modeling Efficiency).
    """
    obs = np.asarray(obs, dtype=float)
    sim = np.asarray(sim, dtype=float)
    
    valid = ~np.isnan(obs) & ~np.isnan(sim)
    obs = obs[valid]
    sim = sim[valid]
    
    n = len(obs)
    if n == 0:
        return np.nan, np.nan, np.nan, 0

    # RMSE
    mse = np.mean((sim - obs) ** 2)
    rmse = np.sqrt(mse)

    # nRMSE (normalisée par la moyenne des observations)
    mean_obs = np.mean(obs)
    nrmse = (rmse / mean_obs) if mean_obs != 0 else np.nan

    # EF (Nash-Sutcliffe / Modeling Efficiency)
    ss_res = np.sum((sim - obs) ** 2)
    ss_tot = np.sum((obs - mean_obs) ** 2)
    ef = (1.0 - (ss_res / ss_tot)) if ss_tot != 0 else np.nan

    return rmse, nrmse, ef, n


def plot_comparison_and_scatter(
    sim_results, bchmk_df, var_name, y_label, title, bchmk_label, 
    output_dir, xlim_bounds, scale_sim=1.0, ylim_Profile=None, legend_loc='upper left',
    metrics_records=None
):
    """
    Génère pour une variable :
      1. Le profil par métamère (Simulations vs Référence)
      2. Le graphique Simulé vs Observé (1:1) avec RMSE, nRMSE et EF dans la légende
    """
    # --------------------------------------------------------------------------
    # GRAPHIQUE 1 : Profil par métamère
    # --------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.set_xlim(xlim_bounds)

    all_y_vals = list(bchmk_df[var_name].dropna() * 1.05)

    for i, (label, df_res) in enumerate(sim_results.items()):
        tmp = df_res[['metamer', var_name]].drop_duplicates().sort_values('metamer').reset_index(drop=True)
        y_vals = tmp[var_name] * scale_sim
        all_y_vals.extend(list(y_vals.dropna() * 1.05))
        
        color = COLORS[i % len(COLORS)]
        marker = MARKERS[i % len(MARKERS)]
        ax.plot(tmp['metamer'], y_vals, color=color, marker=marker, label=label)

    ax.plot(bchmk_df['metamer'], bchmk_df[var_name], color='orange', marker='o', 
            linestyle='--', linewidth=2, label=bchmk_label)

    if ylim_Profile is not None:
        ax.set_ylim(ylim_Profile)
    elif all_y_vals:
        ax.set_ylim(ymin=0, ymax=np.nanmax(all_y_vals))

    ax.set_xlabel('Metamer')
    ax.set_ylabel(y_label)
    ax.set_title(title)
    ax.legend(loc=legend_loc)
    plt.savefig(os.path.join(output_dir, f"{var_name}.PNG"), dpi=200, bbox_inches='tight')
    plt.close()

    # --------------------------------------------------------------------------
    # GRAPHIQUE 2 : Simulé vs Observé (1:1) + Calcul RMSE, nRMSE, EF
    # --------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(6.5, 6.5))
    all_obs_sim_vals = []

    obs_clean = bchmk_df[['metamer', var_name]].dropna().rename(columns={var_name: 'obs'})

    for i, (label, df_res) in enumerate(sim_results.items()):
        sim_clean = df_res[['metamer', var_name]].drop_duplicates().copy()
        sim_clean['sim'] = sim_clean[var_name] * scale_sim
        sim_clean = sim_clean[['metamer', 'sim']].dropna()

        merged = pd.merge(obs_clean, sim_clean, on='metamer', how='inner')
        if merged.empty:
            continue

        rmse, nrmse, ef, n_pts = compute_metrics(merged['obs'], merged['sim'])

        if metrics_records is not None:
            metrics_records.append({
                'Variable': var_name,
                'Reference': bchmk_label,
                'Simulation': label,
                'N': n_pts,
                'RMSE': rmse,
                'nRMSE': nrmse,
                'nRMSE (%)': nrmse * 100 if not np.isnan(nrmse) else np.nan,
                'EF': ef
            })

        all_obs_sim_vals.extend(list(merged['obs']) + list(merged['sim']))

        color = COLORS[i % len(COLORS)]
        marker = MARKERS[i % len(MARKERS)]
        stat_label = f"{label}\n(RMSE={rmse:.2g} | nRMSE={nrmse*100:.1f}% | EF={ef:.2f})"
        
        ax.scatter(merged['obs'], merged['sim'], color=color, marker=marker, s=55, alpha=0.85, label=stat_label)

        for _, row in merged.iterrows():
            ax.annotate(str(int(row['metamer'])), (row['obs'], row['sim']), 
                        fontsize=7, alpha=0.6, xytext=(3, 3), textcoords='offset points')

    if all_obs_sim_vals:
        max_val = np.nanmax(all_obs_sim_vals) * 1.08
        min_val = min(0, np.nanmin(all_obs_sim_vals) * 0.95)
        ax.plot([min_val, max_val], [min_val, max_val], color='black', linestyle='--', linewidth=1.2, label='1:1 line')
        ax.set_xlim(min_val, max_val)
        ax.set_ylim(min_val, max_val)
        ax.set_aspect('equal', adjustable='box')

    ax.set_xlabel(f"Observed - {y_label} ({bchmk_label})")
    ax.set_ylabel(f"Simulated - {y_label}")
    ax.set_title(f"{title} : Simulated vs Observed")
    ax.legend(loc='upper left', fontsize=8.5, framealpha=0.9)
    ax.grid(True, linestyle=':', alpha=0.5)
    plt.savefig(os.path.join(output_dir, f"{var_name}_Sim_vs_Obs.PNG"), dpi=200, bbox_inches='tight')
    plt.close()


def load_and_process_simulation(folder_path, filename=HIDDENZONES_FILENAME):
    """
    Charge et extrait les dimensions finales par métamère pour une simulation donnée.
    """
    filepath = resolve_filepath(folder_path, filename)
    if filepath is None:
        print(f"[Attention] Fichier '{filename}' introuvable dans : {folder_path}")
        return None

    res_raw = pd.read_csv(filepath)
    
    res = res_raw[(res_raw['axis'] == 'MS') & (res_raw['plant'] == 1) & ~np.isnan(res_raw['leaf_Lmax'])].copy()
    if res.empty:
        print(f"[Attention] Aucune donnée valide (MS, plant 1) dans : {filepath}")
        return None

    res_IN = res[~np.isnan(res['internode_Lmax'])].copy()

    last_value_idx = res.groupby(['metamer'])['t'].transform('max') == res['t']
    res = res[last_value_idx].copy()
    res['lamina_Wmax'] = res['leaf_Wmax']
    res['lamina_W_Lg'] = res['leaf_Wmax'] / res['lamina_Lmax']

    if not res_IN.empty:
        last_value_idx_IN = res_IN.groupby(['metamer'])['t'].transform('max') == res_IN['t']
        res_IN = res_IN[last_value_idx_IN].copy()
        res = res[['metamer', 'leaf_Lmax', 'lamina_Lmax', 'sheath_Lmax', 'lamina_Wmax', 'lamina_W_Lg', 'SSLW', 'LSSW']].merge(
            res_IN[['metamer', 'internode_Lmax']], on='metamer', how='outer'
        ).copy()
    else:
        res = res[['metamer', 'leaf_Lmax', 'lamina_Lmax', 'sheath_Lmax', 'lamina_Wmax', 'lamina_W_Lg', 'SSLW', 'LSSW']].copy()
        res['internode_Lmax'] = np.nan

    return res


def load_and_compute_rer(folder_path, hz_filename=HIDDENZONES_FILENAME, axes_filename=AXES_FILENAME):
    """
    Calcule le RER simulé durant la phase quasi-exponentielle pour chaque métamère
    à partir de hiddenzones_outputs.csv et axes_outputs.csv.
    """
    hz_filepath = resolve_filepath(folder_path, hz_filename)
    axes_filepath = resolve_filepath(folder_path, axes_filename)

    if hz_filepath is None or axes_filepath is None:
        print(f"[Attention] Fichiers nécessaires au calcul du RER introuvables dans : {folder_path}")
        return None

    # Import des sorties de simulation
    data_RER = pd.read_csv(hz_filepath)
    if 'plant' in data_RER.columns:
        data_RER = data_RER[(data_RER['axis'] == 'MS') & (data_RER['plant'] == 1) & (data_RER['metamer'] >= 1)].copy()
    else:
        data_RER = data_RER[(data_RER['axis'] == 'MS') & (data_RER['metamer'] >= 1)].copy()
    data_RER.sort_values(['t', 'metamer'], inplace=True)

    data_teq = pd.read_csv(axes_filepath)
    if 'plant' in data_teq.columns:
        data_teq = data_teq[(data_teq['axis'] == 'MS') & (data_teq['plant'] == 1)].copy()
    else:
        data_teq = data_teq[data_teq['axis'] == 'MS'].copy()

    # Temps d'émergence de la feuille précédente
    tmp = data_RER[data_RER['leaf_is_emerged'] == True]
    if tmp.empty:
        return None

    leaf_em = tmp.groupby('metamer', as_index=False)['t'].min()
    leaf_em['t_em'] = leaf_em['t']
    prev_leaf_em = leaf_em.copy()
    prev_leaf_em['metamer'] = leaf_em['metamer'] + 1

    data_RER2 = pd.merge(data_RER, prev_leaf_em[['metamer', 't_em']], on='metamer')
    data_RER2 = data_RER2[data_RER2['t'] <= data_RER2['t_em']]

    # SumTimeEq
    data_teq['SumTimeEq'] = np.cumsum(data_teq['delta_teq'])
    data_RER3 = pd.merge(data_RER2, data_teq[['t', 'SumTimeEq']], on='t')

    # logL
    data_RER3 = data_RER3[data_RER3['leaf_L'] > 0].copy()
    data_RER3['logL'] = np.log(data_RER3['leaf_L'])

    # Estimation du RER par régression OLS
    RER_sim = {}
    for leaf in data_RER3['metamer'].drop_duplicates():
        Y = data_RER3['logL'][data_RER3['metamer'] == leaf]
        X = data_RER3['SumTimeEq'][data_RER3['metamer'] == leaf]
        if len(X) < 2:
            continue
        X_const = sm.add_constant(X)
        mod = sm.OLS(Y, X_const)
        fit_RER = mod.fit()
        if 'SumTimeEq' in fit_RER.params:
            RER_sim[leaf] = fit_RER.params['SumTimeEq']

    if not RER_sim:
        return None

    df_rer = pd.DataFrame(list(RER_sim.items()), columns=['metamer', 'RER']).sort_values('metamer').reset_index(drop=True)
    return df_rer


def plot_rer_comparison_and_scatter(rer_results, data_obs, output_dir, metrics_records=None):
    """
    Trace :
      1. Le graphique RER_comparison.PNG (Simulations + Observations avec barres d'erreur + Paramètres modèle)
      2. Le graphique RER_Sim_vs_Obs.PNG (1:1 avec RMSE, nRMSE et EF)
    """
    # --------------------------------------------------------------------------
    # GRAPHIQUE 1 : RER_comparison.PNG
    # --------------------------------------------------------------------------
    fig, ax1 = plt.subplots(figsize=(7.5, 5))
    ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))

    for i, (label, df_rer) in enumerate(rer_results.items()):
        color = COLORS[i % len(COLORS)]
        marker = MARKERS[i % len(MARKERS)]
        ax1.plot(df_rer['metamer'], df_rer['RER'], label=f"Simulated RER ({label})", 
                 linestyle='-', color=color, marker=marker, markersize=5)

    # Observations Ljutovac 2002 avec barres d'erreur
    if data_obs is not None and 'RER' in data_obs.columns:
        yerr = data_obs['RER_confint'] if 'RER_confint' in data_obs.columns else None
        ax1.errorbar(data_obs['metamer'], data_obs['RER'], yerr=yerr, 
                     marker='o', color='orange', linestyle='', label="Observed RER (Ljutovac 2002)", markersize=5, capsize=3)

    # Paramètres du modèle elongwheat (si le module est installé dans l'environnement)
    if HAS_ELONGWHEAT:
        rer_param = dict((k, v) for k, v in elongwheat_parameters.PARAMETERS.RERmax.items())
        ax1.plot(list(rer_param.keys()), list(rer_param.values()), 
                 marker='*', color='k', linestyle='', markersize=7, label="Model parameters")

    ax1.set_ylabel('Relative Elongation Rate at 12°C (s$^{-1}$)')
    ax1.set_xlabel('Phytomer rank')
    ax1.set_title('RER comparison during exponential-like phase')
    ax1.set_ylim(bottom=0., top=6e-6)
    ax1.set_xlim(left=2)
    ax1.legend(loc='upper left', fontsize=8.5)
    plt.savefig(os.path.join(output_dir, 'RER_comparison.PNG'), format='PNG', bbox_inches='tight', dpi=200)
    plt.close()

    # --------------------------------------------------------------------------
    # GRAPHIQUE 2 : RER_Sim_vs_Obs.PNG + Calcul RMSE, nRMSE, EF
    # --------------------------------------------------------------------------
    if data_obs is None or 'RER' not in data_obs.columns:
        return

    fig, ax = plt.subplots(figsize=(6.5, 6.5))
    ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))

    obs_clean = data_obs[['metamer', 'RER']].dropna().rename(columns={'RER': 'obs'})
    all_obs_sim_vals = []

    for i, (label, df_rer) in enumerate(rer_results.items()):
        sim_clean = df_rer[['metamer', 'RER']].dropna().rename(columns={'RER': 'sim'})
        merged = pd.merge(obs_clean, sim_clean, on='metamer', how='inner')
        if merged.empty:
            continue

        rmse, nrmse, ef, n_pts = compute_metrics(merged['obs'], merged['sim'])

        if metrics_records is not None:
            metrics_records.append({
                'Variable': 'RER',
                'Reference': 'Ljutovac 2002',
                'Simulation': label,
                'N': n_pts,
                'RMSE': rmse,
                'nRMSE': nrmse,
                'nRMSE (%)': nrmse * 100 if not np.isnan(nrmse) else np.nan,
                'EF': ef
            })

        all_obs_sim_vals.extend(list(merged['obs']) + list(merged['sim']))

        color = COLORS[i % len(COLORS)]
        marker = MARKERS[i % len(MARKERS)]
        stat_label = f"{label}\n(RMSE={rmse:.2e} | nRMSE={nrmse*100:.1f}% | EF={ef:.2f})"

        ax.scatter(merged['obs'], merged['sim'], color=color, marker=marker, s=55, alpha=0.85, label=stat_label)

        for _, row in merged.iterrows():
            ax.annotate(str(int(row['metamer'])), (row['obs'], row['sim']), 
                        fontsize=7, alpha=0.6, xytext=(3, 3), textcoords='offset points')

    if all_obs_sim_vals:
        max_val = np.nanmax(all_obs_sim_vals) * 1.08
        min_val = min(0, np.nanmin(all_obs_sim_vals) * 0.95)
        ax.plot([min_val, max_val], [min_val, max_val], color='black', linestyle='--', linewidth=1.2, label='1:1 line')
        ax.set_xlim(min_val, max_val)
        ax.set_ylim(min_val, max_val)
        ax.set_aspect('equal', adjustable='box')

    ax.set_xlabel('Observed - RER at 12°C (s$^{-1}$) (Ljutovac 2002)')
    ax.set_ylabel('Simulated - RER at 12°C (s$^{-1}$)')
    ax.set_title('RER : Simulated vs Observed')
    ax.legend(loc='upper left', fontsize=8.5, framealpha=0.9)
    ax.grid(True, linestyle=':', alpha=0.5)
    plt.savefig(os.path.join(output_dir, 'RER_Sim_vs_Obs.PNG'), format='PNG', bbox_inches='tight', dpi=200)
    plt.close()


# ==============================================================================
# 3. EXECUTION PRINCIPALE
# ==============================================================================

def main():
    plt.ioff()

    sim_results = {}
    rer_results = {}

    for label, folder in SIMULATION_FOLDERS.items():
        df_res = load_and_process_simulation(folder)
        if df_res is not None:
            sim_results[label] = df_res

        df_rer = load_and_compute_rer(folder)
        if df_rer is not None:
            rer_results[label] = df_rer

    if not sim_results and not rer_results:
        print("Erreur : Aucune donnée de simulation n'a pu être chargée. Vérifiez vos chemins.")
        return

    global_min_metamer = min(int(df['metamer'].min()) for df in sim_results.values()) if sim_results else 1
    global_max_metamer = max(int(df['metamer'].max()) for df in sim_results.values()) if sim_results else 12
    xlim_bounds = (global_min_metamer - 1, global_max_metamer + 1)

    metrics_records = []

    # ==========================================================================
    # PARTIE 1 : Comparaison des dimensions avec Ljutovac 2002
    # ==========================================================================
    data_obs = None
    if os.path.isfile(LJUTOVAC_FILEPATH):
        data_obs = pd.read_csv(LJUTOVAC_FILEPATH)
        bchmk_ljutovac = data_obs.loc[data_obs['metamer'] >= global_min_metamer].copy()
        bchmk_ljutovac['lamina_W_Lg'] = bchmk_ljutovac['lamina_Wmax'] / bchmk_ljutovac['lamina_Lmax']

        var_list_cm = ['leaf_Lmax', 'lamina_Lmax', 'sheath_Lmax', 'lamina_Wmax', 'internode_Lmax']
        for var in var_list_cm:
            plot_comparison_and_scatter(
                sim_results=sim_results,
                bchmk_df=bchmk_ljutovac,
                var_name=var,
                y_label=f"{var} (cm)",
                title=var,
                bchmk_label='Ljutovac 2002',
                output_dir=OUTPUT_GRAPHS_DIR,
                xlim_bounds=xlim_bounds,
                scale_sim=100.0,
                metrics_records=metrics_records
            )

        # Ratio lamina_W_Lg
        plot_comparison_and_scatter(
            sim_results=sim_results,
            bchmk_df=bchmk_ljutovac,
            var_name='lamina_W_Lg',
            y_label='lamina_W_Lg',
            title='lamina_W_Lg',
            bchmk_label='Ljutovac 2002',
            output_dir=OUTPUT_GRAPHS_DIR,
            xlim_bounds=xlim_bounds,
            scale_sim=1.0,
            metrics_records=metrics_records
        )
    else:
        print(f"[Attention] Fichier Ljutovac introuvable ({LJUTOVAC_FILEPATH}).")

    # ==========================================================================
    # PARTIE 1bis : Masses structurales vs. adaptation de Bertheloot 2008
    # ==========================================================================

    # A) SSLW Laminae
    bchmk_sslw = pd.DataFrame.from_dict(
        {1: 15, 2: 23, 3: 25, 4: 18, 5: 22, 6: 25, 7: 20, 8: 23, 9: 26, 10: 28, 11: 31}, 
        orient='index'
    ).rename(columns={0: 'SSLW'})
    bchmk_sslw.index.name = 'metamer'
    bchmk_sslw = bchmk_sslw.reset_index()
    bchmk_sslw = bchmk_sslw[bchmk_sslw['metamer'] >= global_min_metamer]

    plot_comparison_and_scatter(
        sim_results=sim_results,
        bchmk_df=bchmk_sslw,
        var_name='SSLW',
        y_label='Structural Specific Lamina Weight (g.m-2)',
        title='Structural Specific Lamina Weight',
        bchmk_label='adapted from Bertheloot 2008',
        output_dir=OUTPUT_GRAPHS_DIR,
        xlim_bounds=xlim_bounds,
        scale_sim=1.0,
        ylim_Profile=(0, 50),
        legend_loc='lower left',
        metrics_records=metrics_records
    )

    # B) LSSW Sheaths
    bchmk_lssw = pd.DataFrame.from_dict(
        {1: 0.08, 2: 0.09, 3: 0.11, 4: 0.18, 5: 0.17, 6: 0.21, 7: 0.24, 8: 0.4, 9: 0.5, 10: 0.55, 11: 0.65}, 
        orient='index'
    ).rename(columns={0: 'LSSW'})
    bchmk_lssw.index.name = 'metamer'
    bchmk_lssw = bchmk_lssw.reset_index()
    bchmk_lssw = bchmk_lssw[bchmk_lssw['metamer'] >= global_min_metamer]

    plot_comparison_and_scatter(
        sim_results=sim_results,
        bchmk_df=bchmk_lssw,
        var_name='LSSW',
        y_label='Lineic Structural Sheath Weight (g.m-1)',
        title='Lineic Structural Sheath Weight',
        bchmk_label='adapted from Bertheloot 2008',
        output_dir=OUTPUT_GRAPHS_DIR,
        xlim_bounds=xlim_bounds,
        scale_sim=1.0,
        ylim_Profile=(0, 0.8),
        legend_loc='upper left',
        metrics_records=metrics_records
    )

    # ==========================================================================
    # PARTIE 3 : RER durant la phase quasi-exponentielle
    # ==========================================================================
    if rer_results:
        plot_rer_comparison_and_scatter(
            rer_results=rer_results,
            data_obs=data_obs,
            output_dir=OUTPUT_GRAPHS_DIR,
            metrics_records=metrics_records
        )

    # ==========================================================================
    # EXPORT ET AFFICHAGE DU RÉCAPITULATIF DES MÉTRIQUES
    # ==========================================================================
    if metrics_records:
        df_metrics = pd.DataFrame(metrics_records)
        csv_path = os.path.join(OUTPUT_GRAPHS_DIR, 'metrics_summary.csv')
        df_metrics.to_csv(csv_path, index=False)
        
        print("\n=== RÉCAPITULATIF DES MÉTRIQUES (RMSE, nRMSE, EF) ===")
        print(df_metrics.to_string(index=False))
        print(f"\nTableau sauvegardé dans : {csv_path}")

    print(f"Terminé ! Tous les graphiques sont dans : '{OUTPUT_GRAPHS_DIR}'")


if __name__ == '__main__':
    main()