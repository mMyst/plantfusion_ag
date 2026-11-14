import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

# --- 1. Paramètres initiaux extraits de la feuille 'En' ---
init_s1 = 0.35    # profilNodeIs1
init_i1 = -0.35   # profilNodeIi1
init_s2 = -0.02   # profilNodeIs2
init_i2 = 1.04    # profilNodeIi2
init_len_max = 5.5 # LenMax
n_noeuds = 20     # Nombre de nœuds (proxy pour le temps)

# Abscisse : rang des nœuds
x = np.arange(1, n_noeuds + 1)

# --- 2. Fonction de calcul basée sur le modèle ---
def calculer_donnees(s1, i1, s2, i2, len_max):
    # Calcul des deux fonctions linéaires
    y1 = s1 * x + i1
    y2 = s2 * x + i2
    
    # Le profil relatif prend le minimum des deux droites, sans descendre sous 0
    profil_relatif = np.maximum(0, np.minimum(y1, y2))
    
    # La hauteur est la somme cumulée de la longueur de chaque entre-nœud
    hauteur_cumulee = np.cumsum(profil_relatif * len_max)
    
    return profil_relatif, hauteur_cumulee

profil_init, hauteur_init = calculer_donnees(init_s1, init_i1, init_s2, init_i2, init_len_max)

# --- 3. Configuration de l'interface graphique ---
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
plt.subplots_adjust(bottom=0.4) # Libère de l'espace en bas pour les curseurs

# Graphique : Profil Relatif
ligne_profil, = ax1.plot(x, profil_init, 'b-o', linewidth=2)
ax1.set_title('Profil Relatif des Nœuds')
ax1.set_xlabel('Rang du nœud')
ax1.set_ylabel('Longueur relative')
ax1.set_ylim(0, 1.5)
ax1.set_xticks(x)
ax1.grid(True, linestyle='--', alpha=0.7)

# Graphique : Hauteur
ligne_hauteur, = ax2.plot(x, hauteur_init, 'r-o', linewidth=2)
ax2.set_title('Hauteur au cours du temps (Cumul)')
ax2.set_xlabel('Rang du nœud (Temps)')
ax2.set_ylabel('Hauteur totale cumulée')
ax2.set_ylim(0, max(hauteur_init) * 1.2 + 2)
ax2.set_xticks(x[::2])
ax2.grid(True, linestyle='--', alpha=0.7)

# --- 4. Création des curseurs (Sliders) ---
couleur_axes = 'lightgray'
ax_s1 = plt.axes([0.15, 0.25, 0.7, 0.03], facecolor=couleur_axes)
ax_i1 = plt.axes([0.15, 0.20, 0.7, 0.03], facecolor=couleur_axes)
ax_s2 = plt.axes([0.15, 0.15, 0.7, 0.03], facecolor=couleur_axes)
ax_i2 = plt.axes([0.15, 0.10, 0.7, 0.03], facecolor=couleur_axes)
ax_len = plt.axes([0.15, 0.05, 0.7, 0.03], facecolor=couleur_axes)

slider_s1 = Slider(ax_s1, 'profilNodeIs1', -1.0, 1.0, valinit=init_s1)
slider_i1 = Slider(ax_i1, 'profilNodeIi1', -2.0, 2.0, valinit=init_i1)
slider_s2 = Slider(ax_s2, 'profilNodeIs2', -1.0, 1.0, valinit=init_s2)
slider_i2 = Slider(ax_i2, 'profilNodeIi2', -2.0, 2.0, valinit=init_i2)
slider_len = Slider(ax_len, 'LenMax', 0.1, 20.0, valinit=init_len_max)

# --- 5. Fonction de mise à jour interactive ---
def actualiser(val):
    profil_rel, hauteur = calculer_donnees(
        slider_s1.val, slider_i1.val, 
        slider_s2.val, slider_i2.val, 
        slider_len.val
    )
    
    # Mise à jour des courbes
    ligne_profil.set_ydata(profil_rel)
    ligne_hauteur.set_ydata(hauteur)
    
    # Ajustement dynamique de l'axe des ordonnées pour la hauteur
    ax2.set_ylim(0, max(hauteur) * 1.1 + 1)
    fig.canvas.draw_idle()

# Connexion des curseurs à la fonction de mise à jour
slider_s1.on_changed(actualiser)
slider_i1.on_changed(actualiser)
slider_s2.on_changed(actualiser)
slider_i2.on_changed(actualiser)
slider_len.on_changed(actualiser)

plt.show()