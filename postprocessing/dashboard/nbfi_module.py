import pandas as pd
import numpy as np
import plotly.graph_objects as go

def _process_experimental_data(ref_path, data_path):
    """Fonction interne pour traiter les données expérimentales."""
    # 1. Chargement du référentiel
    ref_bacs = pd.read_excel(ref_path, sheet_name='ref_bac')
    
    # 2. Chargement de toutes les feuilles de mesures
    all_sheets = pd.read_excel(data_path, sheet_name=None)
    
    # Concaténation des feuilles avec un identifiant
    df_list = []
    for sheet_name, df in all_sheets.items():
        df['sheet_id'] = sheet_name
        df_list.append(df)
    df_final = pd.concat(df_list, ignore_index=True)
    
    # 3. Jointure
    df_final = pd.merge(df_final, ref_bacs, on='ID_BAC', how='left')
    
    # 4. Dates et calcul du DAS
    # Attention au format de la date 'dmy' qui correspond au '%d%m%Y' de R
    df_final['Date'] = pd.to_datetime(df_final['Date'], format='%d%m%Y')
    origin_date = pd.to_datetime("2024-09-30")
    df_final['DAS'] = (df_final['Date'] - origin_date).dt.days
    
    # 5. Filtrage et statistiques (Moyenne et Ecart-type)
    filt_unreg = df_final[df_final['Modality'] == 'LUZ'].copy()
    
    filt_summ = filt_unreg.groupby('DAS').agg(
        NBFI_mean=('NBFI', 'mean'),
        NBFI_sd=('NBFI', 'std')
    ).reset_index()
    
    return filt_summ


def _process_model_data(file_path, model_name):
    """Fonction interne pour traiter un fichier de simulation."""
    # Chargement
    df = pd.read_csv(file_path, sep=';', decimal='.')
    
    # Filtre sur NBI
    nbfi_model = df[df['V1'] == 'NBI'].copy()
    
    # Sélection des colonnes contenant "Timbale"
    timbale_cols = [col for col in nbfi_model.columns if 'Timbale' in col]
    
    # Calculs et création du DataFrame de sortie
    model_processed = pd.DataFrame({
        'steps': nbfi_model['steps'],
        'row_mean': nbfi_model[timbale_cols].mean(axis=1),
        'row_stdev': nbfi_model[timbale_cols].std(axis=1),
        'Model': model_name
    })
    
    # Dates et DAS
    # L'équivalent de as.Date(steps, origin="2024-01-01") en R
    model_processed['Date'] = pd.to_datetime('2024-01-01') + pd.to_timedelta(model_processed['steps'], unit='D')
    model_processed['DAS'] = (model_processed['Date'] - pd.to_datetime("2024-09-30")).dt.days
    
    return model_processed


def NBFI(ref_path, data_path, model_files_dict):
    """
    Fonction principale à appeler dans le dashboard.
    Prend en entrée les chemins des fichiers et retourne un objet graphique Plotly.
    """
    # 1. Traitement des données
    exp_data = _process_experimental_data(ref_path, data_path)
    
    model_dfs = []
    for name, path in model_files_dict.items():
        model_dfs.append(_process_model_data(path, name))
    
    all_models = pd.concat(model_dfs, ignore_index=True)
    
    # 2. Création du graphique interactif avec Plotly
    fig = go.Figure()
    
    # Palette de couleurs pour les modèles
    colors = ['#1f77b4', '#2ca02c', '#ff7f0e', '#9467bd', '#17becf']
    
    # Boucle pour ajouter chaque modèle au graphique
    for i, model_name in enumerate(model_files_dict.keys()):
        m_data = all_models[all_models['Model'] == model_name]
        color = colors[i % len(colors)]
        
        # Ajout du "Ruban" (Intervalle écart-type)
        # Plotly gère cela en traçant la ligne haute, puis la ligne basse à l'envers, et en remplissant l'espace
        fig.add_trace(go.Scatter(
            x=pd.concat([m_data['DAS'], m_data['DAS'][::-1]]),
            y=pd.concat([m_data['row_mean'] + m_data['row_stdev'], 
                         (m_data['row_mean'] - m_data['row_stdev'])[::-1]]),
            fill='toself',
            fillcolor=color,
            opacity=0.2,
            line=dict(color='rgba(255,255,255,0)'), # Ligne invisible
            hoverinfo="skip",
            showlegend=False,
            name=f"{model_name} SD"
        ))
        
        # Ajout de la ligne moyenne
        fig.add_trace(go.Scatter(
            x=m_data['DAS'],
            y=m_data['row_mean'],
            mode='lines',
            line=dict(color=color, width=2),
            name=model_name
        ))

    # 3. Ajout des données expérimentales observées
    fig.add_trace(go.Scatter(
        x=exp_data['DAS'],
        y=exp_data['NBFI_mean'],
        mode='lines+markers',
        line=dict(color='red', width=2, dash='dash'),
        marker=dict(color='red', size=8),
        error_y=dict(
            type='data',
            array=exp_data['NBFI_sd'],
            visible=True,
            color='red',
            thickness=1.5
        ),
        name='Observations (LUZ)'
    ))
    
    # 4. Mise en page du graphique
    fig.update_layout(
        title='Observed vs Simulated Height (NBFI)',
        xaxis_title='Days After Sowing (DAS)',
        yaxis_title='NBFI',
        template='plotly_white', # Fond blanc propre
        hovermode="x unified",   # Affiche toutes les valeurs au survol d'un DAS
        legend=dict(
            orientation="h",     # Légende horizontale sous le graph
            yanchor="bottom",
            y=-0.3,
            xanchor="center",
            x=0.5
        )
    )
    
    return fig