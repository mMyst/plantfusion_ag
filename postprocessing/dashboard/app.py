import streamlit as st
from nbfi_module import NBFI

# Configuration de la page (doit être la première commande Streamlit)
st.set_page_config(
    page_title="Dashboard LUBBAC",
    page_icon="🌱",
    layout="wide"
)

# --- TITRE PRINCIPAL ---
st.title("🌱 Comparaison des Modèles - Projet LUBBAC")
st.markdown("Ce tableau de bord permet de comparer les données observées sur la luzerne avec les différentes sorties de simulation.")
st.divider()

# --- BARRE LATÉRALE (SIDEBAR) ---
st.sidebar.header("Paramètres")

# 1. Chemins des données expérimentales (fixes pour l'instant)
ref_path = r"C:\Users\agrumel\Documents\Données\LUBBAC\ref_LUBBAC.xlsx"
data_path = r"\\pnas1.stockage.inra.fr\nap-urp3f\Projets\2024_Lubbac_GLouarn\Mesures_Morpho_Luzerne LUBBAC_avecLignesGrises.xlsx"

# 2. Dictionnaire de tous les modèles disponibles
tous_les_modeles = {
    "Modèle 1 (Standard)": r"C:\Users\agrumel\Desktop\toto_2_l-egume_Timbale_LUBBAC-Timbale_LUBBAC_random88_scenario1-1_LUBBAC_nomgmt_0_LUBBAC_D_24_25_SD9-9_.csv",
    "Modèle 2 (Couplage)": r"C:\Users\agrumel\Documents\Données\Sorties Couplage\full_coupling_LUBBAC_1900_nogel\legume\brut\2_Timbale_LUBBAC__LUBBAC_nomgmt_0_LUBBAC_D_24_25_SD9-9_\toto_2_l-egume_Timbale_LUBBAC-Timbale_LUBBAC_random88_scenario1-1_LUBBAC_nomgmt_0_LUBBAC_D_24_25_SD9-9_.csv",
    "Modèle 3 (L-egume)": r"C:\Users\agrumel\Documents\Données\Sorties L-egume\legume_LUBBAC\legume\brut\legume_lubbac\toto_2_l-egume_Timbale_LUBBAC-Timbale_LUBBAC_random88_scenario1-1_LUBBAC_nomgmt_0_LUBBAC_D_24_25_SD9-9_.csv"
}

# 3. Widget interactif : choix des modèles à afficher
modeles_selectionnes = st.sidebar.multiselect(
    "Quels modèles voulez-vous afficher ?",
    options=list(tous_les_modeles.keys()),
    default=list(tous_les_modeles.keys()) # Par défaut, on coche tout
)

# On filtre notre dictionnaire en fonction de ce que l'utilisateur a coché
modeles_a_comparer = {nom: chemin for nom, chemin in tous_les_modeles.items() if nom in modeles_selectionnes}

# --- AFFICHAGE DU GRAPHIQUE ---
st.subheader("Évolution du NBFI (Nombre de Feuilles)")

# On vérifie qu'au moins un modèle est sélectionné
if len(modeles_a_comparer) > 0:
    # Optionnel : ajouter un petit message de chargement
    with st.spinner('Génération du graphique en cours...'):
        # On appelle ta fonction Python !
        fig_nbfi = NBFI(ref_path, data_path, modeles_a_comparer)
        
        # On affiche la figure Plotly dans Streamlit
        st.plotly_chart(fig_nbfi, use_container_width=True)
else:
    st.warning("👈 Veuillez sélectionner au moins un modèle dans le menu de gauche.")