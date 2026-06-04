#this is a script that will be used to create a dashboard for the outputs of LUBBAC simulations

#we will build interactible graphs using plotly, and a recap pdf of the main graphs

#there's going to be 3 main sections : 
# - outputs from wheat 
# - outputs from legume
# - outputs from soil

#we will compare simulated data with measurements from the LUBBAC trial, such as : 
# plant height, leaf apparition, plant N content, intercepted PAR, transmitted PAR,  soil N

#the reference data will always be the same, what can change is the output files to compare with



###path for wheat reference data
#height
wheat_ref_height = r"\\pnas1.stockage.inra.fr\nap-urp3f\Projets\2024_Lubbac_GLouarn\Developpement_ble.xlsx",sheet = Hauteurs
#leaf apparition
wheat_ref_phyllochron = r"\\pnas1.stockage.inra.fr\nap-urp3f\Projets\2024_Lubbac_GLouarn\Developpement_ble.xlsx",sheet = data
#N content

###path for legume reference data
#height
#legume_ref_height =

#leaf apparition
#N content

###path for environment reference data
#transmitted PAR


# Fichier : mon_dashboard.py

# On importe la fonction depuis l'autre fichier
from nbfi_module import NBFI 

# 1. Définir les chemins des données d'observation
ref_path = r"C:\Users\agrumel\Documents\Données\LUBBAC\ref_LUBBAC.xlsx"
data_path = r"\\pnas1.stockage.inra.fr\nap-urp3f\Projets\2024_Lubbac_GLouarn\Mesures_Morpho_Luzerne LUBBAC_avecLignesGrises.xlsx"

# 2. Définir un dictionnaire avec le Nom du modèle et son chemin CSV
# Tu peux en ajouter autant que tu veux sans modifier la fonction NBFI() !
modeles_a_comparer = {
    "Modèle 1 (Standard)": r"C:\Users\agrumel\Desktop\toto_2_l-egume_Timbale_LUBBAC..._.csv",
    "Modèle 2 (Couplage)": r"C:\Users\agrumel\Documents\Données\Sorties Couplage\...\toto_2..._.csv",
    "Modèle 3 (L-egume)": r"C:\Users\agrumel\Documents\Données\Sorties L-egume\...\toto_2..._.csv"
}

# 3. Générer la figure
fig_nbfi = NBFI(ref_path, data_path, modeles_a_comparer)

# 4. Afficher la figure (Dépend du framework de ton dashboard)
# Si tu testes juste dans un notebook ou script simple : 
# fig_nbfi.show()

# Si tu utilises Streamlit :
# import streamlit as st
# st.plotly_chart(fig_nbfi, use_container_width=True)

# Si tu utilises Dash :
# from dash import dcc
# dcc.Graph(figure=fig_nbfi)