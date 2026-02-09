import os
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

delta_t_simuls = 0 #1509
meteo_data = pd.read_csv(os.path.join('inputs_fspmwheat', 'meteo_Ljutovac2002.csv'), index_col='t')
meteo_data['Date'] = pd.to_datetime(meteo_data['Date'], format='%d/%m/%Y')


def phloem(df_current_organs, df_marion_organs, tmin, tmax):
    fig, axs = plt.subplots(2, 2)

    # phloem sucrose & AA

    # 1
    axs[0, 0].plot(df_current_organs[(df_current_organs.organ == 'phloem')]['t'], df_current_organs[(df_current_organs.organ == 'phloem')]['Conc_Sucrose'], label='current')
    axs[0, 0].plot(df_marion_organs[(df_marion_organs.organ == 'phloem')]['t'], df_marion_organs[(df_marion_organs.organ == 'phloem')]['Conc_Sucrose'], label='Control')
    axs[0, 0].legend()
    axs[0, 0].set_xlim(tmin, tmax)
    axs[0, 0].set_ylabel('Concentration sucrose (µmol g-1)')

    # 2
    axs[0, 1].plot(df_current_organs[(df_current_organs.organ == 'phloem')]['t'], df_current_organs[(df_current_organs.organ == 'phloem')]['sucrose'], label='current')
    axs[0, 1].plot(df_marion_organs[(df_marion_organs.organ == 'phloem')]['t'], df_marion_organs[(df_marion_organs.organ == 'phloem')]['sucrose'], label='Control')
    axs[0, 1].legend()
    axs[0, 1].set_xlim(tmin, tmax)
    # axs[0, 1].set_ylim(0, 500)
    axs[0, 1].set_ylabel('Amount of sucrose (µmol C)')

    # 3
    axs[1, 0].plot(df_current_organs[(df_current_organs.organ == 'phloem')]['t'], df_current_organs[(df_current_organs.organ == 'phloem')]['Conc_Amino_Acids'], label='current')
    axs[1, 0].plot(df_marion_organs[(df_marion_organs.organ == 'phloem')]['t'], df_marion_organs[(df_marion_organs.organ == 'phloem')]['Conc_Amino_Acids'], label='Control')
    axs[1, 0].legend()
    axs[1, 0].set_xlim(tmin, tmax)
    # axs[1, 0].set_ylim(0, 200)
    axs[1, 0].set_ylabel('Concentration amino acids (µmol g-1)')

    ax2 = axs[1, 0].twiny()
    ax2.set_xticks(axs[1, 0].get_xticks())
    ax2.set_xticklabels(meteo_data.loc[axs[1, 0].get_xticks()]['Date'].dt.strftime('%d/%m'))
    ax2.xaxis.set_ticks_position('bottom')  # set the position of the second x-axis to bottom
    ax2.xaxis.set_label_position('bottom')  # set the position of the second x-axis to bottom
    ax2.spines['bottom'].set_position(('outward', 35))

    # 4
    axs[1, 1].plot(df_current_organs[(df_current_organs.organ == 'phloem')]['t'], df_current_organs[(df_current_organs.organ == 'phloem')]['amino_acids'], label='current')
    axs[1, 1].plot(df_marion_organs[(df_marion_organs.organ == 'phloem')]['t'], df_marion_organs[(df_marion_organs.organ == 'phloem')]['amino_acids'], label='Control')
    axs[1, 1].legend()
    axs[1, 1].set_xlim(tmin, tmax)
    # axs[1, 1].set_ylim(0, 25)
    axs[1, 1].set_ylabel('amino acids (µmol N)')

    ax2 = axs[1, 1].twiny()
    ax2.set_xticks(axs[1, 1].get_xticks())
    ax2.set_xticklabels(meteo_data.loc[axs[1, 1].get_xticks()]['Date'].dt.strftime('%d/%m'))
    ax2.xaxis.set_ticks_position('bottom')  # set the position of the second x-axis to bottom
    ax2.xaxis.set_label_position('bottom')  # set the position of the second x-axis to bottom
    ax2.spines['bottom'].set_position(('outward', 35))

    plt.tight_layout()
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()

def photosynthesis(df_current_axes, df_marion_axes):
    df_current_axes['day'] = df_current_axes['t'] // 24 +1
    df_marion_axes['day'] = df_marion_axes['t'] // 24 +1

    fig, axis = plt.subplots()
    axis.plot(df_current_axes['day'].unique(), df_current_axes.groupby('day')['Total_Photosynthesis'].sum(), label='current')
    axis.plot(df_marion_axes['day'].unique(), df_marion_axes.groupby('day')['Total_Photosynthesis'].sum(), label='Control')

    axis.set_xlabel('Time (day)')
    axis.set_ylabel('Total Photosynthesis µmol C')
    axis.legend()

    plt.tight_layout()
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()

def roots(df_current_organs, df_marion_organs, tmin, tmax):
    fig, axs = plt.subplots(2, 2)

    # phloem sucrose & AA

    # 1
    axs[0, 0].plot(df_current_organs[(df_current_organs.organ == 'roots')]['t'], df_current_organs[(df_current_organs.organ == 'roots')]['Conc_Sucrose'], label='current')
    axs[0, 0].plot(df_marion_organs[(df_marion_organs.organ == 'roots')]['t'], df_marion_organs[(df_marion_organs.organ == 'roots')]['Conc_Sucrose'], label='Control')
    axs[0, 0].legend()
    axs[0, 0].set_xlim(tmin, tmax)
    axs[0, 0].set_ylabel('Concentration sucrose (µmol g-1)')

    # 2
    axs[0, 1].plot(df_current_organs[(df_current_organs.organ == 'roots')]['t'], df_current_organs[(df_current_organs.organ == 'roots')]['sucrose'], label='current')
    axs[0, 1].plot(df_marion_organs[(df_marion_organs.organ == 'roots')]['t'], df_marion_organs[(df_marion_organs.organ == 'roots')]['sucrose'], label='Control')
    axs[0, 1].legend()
    axs[0, 1].set_xlim(tmin, tmax)
    axs[0, 1].set_ylim(0, 1000)
    axs[0, 1].set_ylabel('Amount of sucrose (µmol C)')

    # 3
    axs[1, 0].plot(df_current_organs[(df_current_organs.organ == 'roots')]['t'], df_current_organs[(df_current_organs.organ == 'roots')]['Conc_Nitrates'], label='current')
    axs[1, 0].plot(df_marion_organs[(df_marion_organs.organ == 'roots')]['t'], df_marion_organs[(df_marion_organs.organ == 'roots')]['Conc_Nitrates'], label='Control')
    axs[1, 0].legend()
    axs[1, 0].set_xlim(tmin, tmax)
    # axs[1, 0].set_ylim(0, 200)
    axs[1, 0].set_ylabel('Concentration nitrates (µmol g-1)')

    ax2 = axs[1, 0].twiny()
    ax2.set_xticks(axs[1, 0].get_xticks())
    ax2.set_xticklabels(meteo_data.loc[axs[1, 0].get_xticks()]['Date'].dt.strftime('%d/%m'))
    ax2.xaxis.set_ticks_position('bottom')  # set the position of the second x-axis to bottom
    ax2.xaxis.set_label_position('bottom')  # set the position of the second x-axis to bottom
    ax2.spines['bottom'].set_position(('outward', 35))

    # 4
    axs[1, 1].plot(df_current_organs[(df_current_organs.organ == 'roots')]['t'], df_current_organs[(df_current_organs.organ == 'roots')]['Conc_cytokinins'], label='current')
    axs[1, 1].plot(df_marion_organs[(df_marion_organs.organ == 'roots')]['t'], df_marion_organs[(df_marion_organs.organ == 'roots')]['Conc_cytokinins'], label='Control')
    axs[1, 1].legend()
    axs[1, 1].set_xlim(tmin, tmax)
    axs[1, 1].set_ylabel('Conc_cytokinins (AU g-1)')

    ax2 = axs[1, 1].twiny()
    ax2.set_xticks(axs[1, 1].get_xticks())
    ax2.set_xticklabels(meteo_data.loc[axs[1, 1].get_xticks()]['Date'].dt.strftime('%d/%m'))
    ax2.xaxis.set_ticks_position('bottom')  # set the position of the second x-axis to bottom
    ax2.xaxis.set_label_position('bottom')  # set the position of the second x-axis to bottom
    ax2.spines['bottom'].set_position(('outward', 35))

    plt.tight_layout()
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()

def dry_mass(df_current_axes, df_marion_axes, df_current_organs, df_marion_organs, tmin, tmax):
    fig, axs = plt.subplots(2, 2, sharex=True)

    # Dry mass shoot
    axs[0, 0].plot(df_current_axes['t'], df_current_axes['sum_dry_mass_shoot'], label='current')
    axs[0, 0].plot(df_marion_axes['t'], df_marion_axes['sum_dry_mass_shoot'], label='Control')
    axs[0, 0].legend()
    axs[0, 0].set_xlim(tmin, tmax)
    # axs[0, 0].set_ylim(0, 1)
    axs[0, 0].set_ylabel('Dry mass shoot (g)')

    # Dry mass roots
    axs[0, 1].plot(df_current_axes['t'], df_current_axes['sum_dry_mass_roots'], label='current')
    axs[0, 1].plot(df_marion_axes['t'], df_marion_axes['sum_dry_mass_roots'], label='Control')
    axs[0, 1].legend()
    axs[0, 1].set_xlim(tmin, tmax)
    # axs[0, 1].set_ylim(0, 1)
    axs[0, 1].set_ylabel('Dry mass roots (g)')

    # mstruct shoot
    axs[1, 0].plot(df_current_axes['t'], df_current_axes['mstruct_shoot'], label='current')
    axs[1, 0].plot(df_marion_axes['t'], df_marion_axes['mstruct_shoot'], label='Control')
    axs[1, 0].legend()
    axs[1, 0].set_xlim(tmin, tmax)
    # axs[1, 0].set_ylim(0, 1)
    axs[1, 0].set_ylabel('mstruct shoot (g)')
    ax2 = axs[1, 0].twiny()
    ax2.set_xticks(axs[1, 0].get_xticks())
    ax2.set_xticklabels(meteo_data.loc[axs[1, 0].get_xticks()]['Date'].dt.strftime('%d/%m'))
    ax2.xaxis.set_ticks_position('bottom')  # set the position of the second x-axis to bottom
    ax2.xaxis.set_label_position('bottom')  # set the position of the second x-axis to bottom
    ax2.spines['bottom'].set_position(('outward', 35))

    # mstruct roots
    axs[1, 1].plot(df_current_organs[df_current_organs['organ'] == 'roots']['t'], df_current_organs[df_current_organs['organ'] == 'roots']['mstruct'], label='current')
    axs[1, 1].plot(df_marion_organs[df_marion_organs['organ'] == 'roots']['t'], df_marion_organs[df_marion_organs['organ'] == 'roots']['mstruct'], label='Control')
    axs[1, 1].legend()
    axs[1, 1].set_xlim(tmin, tmax)
    # axs[1, 1].set_ylim(0, 0.75)
    axs[1, 1].set_ylabel('mstruct roots (g)')
    ax2 = axs[1, 1].twiny()
    ax2.set_xticks(axs[1, 1].get_xticks())
    ax2.set_xticklabels(meteo_data.loc[axs[1, 1].get_xticks()]['Date'].dt.strftime('%d/%m'))
    ax2.xaxis.set_ticks_position('bottom')  # set the position of the second x-axis to bottom
    ax2.xaxis.set_label_position('bottom')  # set the position of the second x-axis to bottom
    ax2.spines['bottom'].set_position(('outward', 35))

    plt.tight_layout()
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()

    # shoot : root
    fig, axis = plt.subplots()
    axis.plot(df_current_axes['t'], df_current_axes['shoot_roots_ratio'], label='current')
    axis.plot(df_marion_axes['t'], df_marion_axes['shoot_roots_ratio'], label='Control')
    axis.legend()
    axis.set_xlim(tmin, tmax)
    axis.set_ylim(0, 2)
    axis.set_ylabel('shoot : root ratio')

    ax2 = axis.twiny()
    ax2.set_xticks(axis.get_xticks())
    ax2.set_xticklabels(meteo_data.loc[axis.get_xticks()]['Date'].dt.strftime('%d/%m'))
    ax2.xaxis.set_ticks_position('bottom')  # set the position of the second x-axis to bottom
    ax2.xaxis.set_label_position('bottom')  # set the position of the second x-axis to bottom
    ax2.spines['bottom'].set_position(('outward', 35))
    plt.tight_layout()
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()


def N_mass(df_current_axes, df_marion_axes, df_current_organs, df_marion_organs, tmin, tmax):
    fig, axs = plt.subplots(2, 2, sharex=True)

    # % N axis
    axs[0, 0].plot(df_current_axes['t'], df_current_axes['N_content'], label='current')
    axs[0, 0].plot(df_marion_axes['t'], df_marion_axes['N_content'], label='Control')
    axs[0, 0].legend()
    axs[0, 0].set_xlim(tmin, tmax)
    axs[0, 0].set_ylim(0, 10)
    axs[0, 0].set_ylabel('N content axis (% DM)')

    # N shoot
    axs[0, 1].plot(df_current_axes['t'], df_current_axes['N_content_shoot'], label='current')
    axs[0, 1].plot(df_marion_axes['t'], df_marion_axes['N_content_shoot'], label='Control')
    axs[0, 1].legend()
    axs[0, 1].set_xlim(tmin, tmax)
    axs[0, 1].set_ylabel('N content shoot (% DM)')

    # N axis
    axs[1, 0].plot(df_current_axes['t'], df_current_axes['sum_N_g'], label='current')
    axs[1, 0].plot(df_marion_axes['t'], df_marion_axes['sum_N_g'], label='Control')
    axs[1, 0].legend()
    axs[1, 0].set_xlim(tmin, tmax)
    axs[1, 0].set_ylabel('N content axis (g)')
    ax2 = axs[1, 0].twiny()
    ax2.set_xticks(axs[1, 0].get_xticks())
    ax2.set_xticklabels(meteo_data.loc[axs[1, 0].get_xticks()]['Date'].dt.strftime('%d/%m'))
    ax2.xaxis.set_ticks_position('bottom')  # set the position of the second x-axis to bottom
    ax2.xaxis.set_label_position('bottom')  # set the position of the second x-axis to bottom
    ax2.spines['bottom'].set_position(('outward', 35))

    # N uptake
    axs[1, 1].plot(df_current_organs[df_current_organs['organ'] == 'roots']['t'], df_current_organs[df_current_organs['organ'] == 'roots']['Uptake_Nitrates'], label='current')
    axs[1, 1].plot(df_marion_organs[df_marion_organs['organ'] == 'roots']['t'], df_marion_organs[df_marion_organs['organ'] == 'roots']['Uptake_Nitrates'], label='Control')
    axs[1, 1].legend()
    axs[1, 1].set_xlim(tmin, tmax)
    axs[1, 1].set_ylabel('Nitrate uptake (µmol)')
    ax2 = axs[1, 1].twiny()
    ax2.set_xticks(axs[1, 1].get_xticks())
    ax2.set_xticklabels(meteo_data.loc[axs[1, 1].get_xticks()]['Date'].dt.strftime('%d/%m'))
    ax2.xaxis.set_ticks_position('bottom')  # set the position of the second x-axis to bottom
    ax2.xaxis.set_label_position('bottom')  # set the position of the second x-axis to bottom
    ax2.spines['bottom'].set_position(('outward', 35))

    plt.tight_layout()
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()


def surface(df_current_elements, df_marion_elements, tmin, tmax):
    fig, axs = plt.subplots(2, 2)

    # Total green area
    axs[0, 0].plot(df_current_elements['t'].unique(), df_current_elements.groupby('t')['green_area'].sum(), label='current')
    axs[0, 0].plot(df_marion_elements['t'].unique(), df_marion_elements.groupby('t')['green_area'].sum(), label='Control')
    axs[0, 0].legend()
    axs[0, 0].set_xlim(tmin, tmax)
    axs[0, 0].set_ylabel('Total green area (m²)')

    # Blade green area
    df_current_elements_blade = df_current_elements[df_current_elements.organ == 'blade']
    df_marion_elements_blade = df_marion_elements[df_marion_elements.organ == 'blade']
    axs[0, 1].plot(df_current_elements_blade['t'].unique(), df_current_elements_blade.groupby('t')['green_area'].sum(), label='current')
    axs[0, 1].plot(df_marion_elements_blade['t'].unique(), df_marion_elements_blade.groupby('t')['green_area'].sum(), label='Control')
    axs[0, 1].legend()
    axs[0, 1].set_xlim(tmin, tmax)
    axs[0, 1].set_ylabel('Blade green area (m²)')

    # Sheath green area
    df_current_elements_sheath = df_current_elements[df_current_elements.organ == 'sheath']
    df_marion_elements_sheath = df_marion_elements[df_marion_elements.organ == 'sheath']
    axs[1, 0].plot(df_current_elements_sheath['t'].unique(), df_current_elements_sheath.groupby('t')['green_area'].sum(), label='current')
    axs[1, 0].plot(df_marion_elements_sheath['t'].unique(), df_marion_elements_sheath.groupby('t')['green_area'].sum(), label='Control')
    axs[1, 0].legend()
    axs[1, 0].set_xlim(tmin, tmax)
    axs[1, 0].set_ylabel('Sheath green area (m²)')
    ax2 = axs[1, 0].twiny()
    ax2.set_xticks(axs[1, 0].get_xticks())
    ax2.set_xticklabels(meteo_data.loc[axs[1, 0].get_xticks()]['Date'].dt.strftime('%d/%m'))
    ax2.xaxis.set_ticks_position('bottom')  # set the position of the second x-axis to bottom
    ax2.xaxis.set_label_position('bottom')  # set the position of the second x-axis to bottom
    ax2.spines['bottom'].set_position(('outward', 35))


    # Internode green area
    df_current_elements_internode = df_current_elements[df_current_elements.organ == 'internode']
    df_marion_elements_internode = df_marion_elements[df_marion_elements.organ == 'internode']
    axs[1, 1].plot(df_current_elements_internode['t'].unique(), df_current_elements_internode.groupby('t')['green_area'].sum(), label='current')
    axs[1, 1].plot(df_marion_elements_internode['t'].unique(), df_marion_elements_internode.groupby('t')['green_area'].sum(), label='Control')
    axs[1, 1].legend()
    axs[1, 1].set_xlim(tmin, tmax)
    axs[1, 1].set_ylabel('Internode green area (m²)')
    ax2 = axs[1, 1].twiny()
    ax2.set_xticks(axs[1, 1].get_xticks())
    ax2.set_xticklabels(meteo_data.loc[axs[1, 1].get_xticks()]['Date'].dt.strftime('%d/%m'))
    ax2.xaxis.set_ticks_position('bottom')  # set the position of the second x-axis to bottom
    ax2.xaxis.set_label_position('bottom')  # set the position of the second x-axis to bottom
    ax2.spines['bottom'].set_position(('outward', 35))

    plt.tight_layout()
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()


def include_images(graphs_path, graphs_Marion_path):
    graphs_to_include = ['leaf_L_hz.PNG', 'Leaf_Lmax.PNG', 'RER_comparison.PNG', 'phyllochron.PNG', 'lamina_Wmax.PNG', 'SSLW.PNG']
    for graph in graphs_to_include:
        im = plt.imread(os.path.join(graphs_path, graph))
        fig = plt.figure(figsize=(13, 10))
        fig.figimage(im)
        pdf.savefig()  # saves the current figure into a pdf page
        plt.close()


def leaf_length_mstruct_area(df_current_hz, df_marion_hz, df_current_elements, df_marion_elements, tmin, tmax):
    # Loop through phytomer
    for phyto_id in df_current_hz['metamer'].unique():
        fig, axs = plt.subplots(3, 3)

        # Length leaf
        axs[0][0].plot(df_current_hz[(df_current_hz.metamer == phyto_id)]['t'], df_current_hz[(df_current_hz.metamer == phyto_id)]['leaf_L'], color='c', label='current')
        if phyto_id in (1,2):
            axs[0][0].plot(df_marion_elements[(df_marion_elements.metamer == phyto_id)]['t'].unique(), df_marion_elements[(df_marion_elements.metamer == phyto_id)].groupby('t')['length'].sum(), color='orange', label='Control')
        else:
            axs[0][0].plot(df_marion_hz[(df_marion_hz.metamer == phyto_id)]['t'], df_marion_hz[(df_marion_hz.metamer == phyto_id)]['leaf_L'], color='orange', label='Control')
        axs[0][0].set_xlim(tmin, tmax)
        axs[0][0].set_ylabel('Leaf_L  hz' + ' (m)')
        axs[0][0].set_title('Leaf' + '_' + str(phyto_id))
        axs[0][0].set_xticks([])

        # Length sheath
        axs[0][1].plot(df_current_elements[(df_current_elements.metamer == phyto_id) & (df_current_elements.organ == 'sheath')]['t'].unique(), df_current_elements[(df_current_elements.metamer == phyto_id) & (df_current_elements.organ == 'sheath')].groupby('t')['length'].sum(), color='c', label='current')
        axs[0][1].plot(df_marion_elements[(df_marion_elements.metamer == phyto_id) & (df_marion_elements.organ == 'sheath')]['t'].unique(), df_marion_elements[(df_marion_elements.metamer == phyto_id) & (df_marion_elements.organ == 'sheath')].groupby('t')['length'].sum(), color='orange', label='Control')
        axs[0][1].set_xlim(tmin, tmax)
        axs[0][1].set_ylabel('Sheath_L elt' + ' (m)')
        axs[0][1].set_xticks([])

        # Mstruct hz
        if phyto_id in (1, 2):
            axs[1][0].plot(df_current_elements[(df_current_elements.metamer == phyto_id) & (df_current_elements.element == 'LeafElement1')]['t'].unique(), df_current_elements[(df_current_elements.metamer == phyto_id) & (df_current_elements.element == 'LeafElement1')]['mstruct'], color='c', label='current')
            axs[1][0].plot(df_marion_elements[(df_marion_elements.metamer == phyto_id)]['t'].unique(), df_marion_elements[(df_marion_elements.metamer == phyto_id) & (df_marion_elements.element == 'LeafElement1')]['mstruct'],
                        color='orange', label='Control')
        else:
            axs[1][0].plot(df_current_hz[(df_current_hz.metamer == phyto_id)]['t'], df_current_hz[(df_current_hz.metamer == phyto_id)]['mstruct'], color='c', label='current')
            axs[1][0].plot(df_marion_hz[(df_marion_hz.metamer == phyto_id)]['t'], df_marion_hz[(df_marion_hz.metamer == phyto_id)]['mstruct'], color='orange', label='Control')
        axs[1][0].set_xlim(tmin, tmax)
        axs[1][0].set_ylabel('mstruct hz' + ' (g)')
        axs[1][0].set_xticks([])

        # mstruct lamina
        if phyto_id == 0: continue
        axs[1][1].plot(df_current_elements[(df_current_elements.metamer == phyto_id) & (df_current_elements.organ == 'blade')]['t'].unique(), df_current_elements[(df_current_elements.metamer == phyto_id) & (df_current_elements.organ == 'blade')].groupby('t')['mstruct'].sum(), color='c', label='current')
        axs[1][1].plot(df_marion_elements[(df_marion_elements.metamer == phyto_id) & (df_marion_elements.organ == 'blade')]['t'].unique(), df_marion_elements[(df_marion_elements.metamer == phyto_id) & (df_marion_elements.organ == 'blade')].groupby('t')['mstruct'].sum(), color='orange', label='Control')
        axs[1][1].set_xlim(tmin, tmax)
        axs[1][1].set_ylabel('Lam mstruct elt' + ' (g)')
        axs[1][1].set_xticks([])

        # mstruct sheath
        axs[2][0].plot(df_current_elements[(df_current_elements.metamer == phyto_id) & (df_current_elements.organ == 'sheath')]['t'].unique(), df_current_elements[(df_current_elements.metamer == phyto_id) & (df_current_elements.organ == 'sheath')].groupby('t')['mstruct'].sum(), color='c', label='current')
        axs[2][0].plot(df_marion_elements[(df_marion_elements.metamer == phyto_id) & (df_marion_elements.organ == 'sheath')]['t'].unique(), df_marion_elements[(df_marion_elements.metamer == phyto_id) & (df_marion_elements.organ == 'sheath')].groupby('t')['mstruct'].sum(), color='orange', label='Control')
        axs[2][0].set_xlim(tmin, tmax)
        axs[2][0].set_ylabel('Sheath mstruct' + ' (g)')

        # Green area lamina
        current_blade_green_area = df_current_elements[(df_current_elements.metamer == phyto_id) & (df_current_elements.element == 'LeafElement1')]
        Marion_blade_green_area = df_marion_elements[(df_marion_elements.metamer == phyto_id) & (df_marion_elements.element == 'LeafElement1')]
        axs[2][1].plot(current_blade_green_area['t'], current_blade_green_area['green_area'], color='c', label='current')
        axs[2][1].plot(Marion_blade_green_area['t'], Marion_blade_green_area['green_area'], color='orange', label='Control')
        axs[2][1].set_xlim(tmin, tmax)
        axs[2][1].set_ylabel('Lamina GA' + ' (m2)')

        # Internode length
        current_internode_L = df_current_hz[df_current_hz.metamer == phyto_id]['internode_L']
        Marion_internode_L = df_marion_hz[df_marion_hz.metamer == phyto_id]['internode_L']
        axs[0][2].plot(df_current_hz[(df_current_hz.metamer == phyto_id)]['t'], current_internode_L, color='c', label='current')
        axs[0][2].plot(df_marion_hz[(df_marion_hz.metamer == phyto_id)]['t'], Marion_internode_L, color='orange', label='Control')
        axs[0][2].set_xlim(tmin, tmax)
        axs[0][2].set_ylabel('Internode_L hz (m)')
        axs[0][2].legend(loc='upper center', bbox_to_anchor=(0.25, 1.5), ncol=2, fontsize="8")

        # N content lamina
        if phyto_id == 0: continue
        current_lamina_L = df_current_elements[(df_current_elements.metamer == phyto_id) & (df_current_elements.element == 'LeafElement1')]['N_tot']
        Marion_lamina_L = df_marion_elements[(df_marion_elements.metamer == phyto_id) & (df_marion_elements.element == 'LeafElement1')]['N_tot']
        axs[1][2].plot(df_current_elements[(df_current_elements.metamer == phyto_id) & (df_current_elements.element == 'LeafElement1')]['t'].unique(), current_lamina_L, color='c', label='current')
        axs[1][2].plot(df_marion_elements[(df_marion_elements.metamer == phyto_id) & (df_marion_elements.element == 'LeafElement1')]['t'].unique(), Marion_lamina_L, color='orange', label='Control')
        axs[1][2].set_xlim(tmin, tmax)
        axs[1][2].set_ylabel('N content lam (g)')

        # N content sheath
        if phyto_id == 0: continue
        axs[2][2].plot(df_current_elements[(df_current_elements.metamer == phyto_id) & (df_current_elements.organ == 'sheath')]['t'].unique(), df_current_elements[(df_current_elements.metamer == phyto_id) & (df_current_elements.organ == 'sheath')].groupby('t')['N_tot'].sum(), color='c', label='current')
        axs[2][2].plot(df_marion_elements[(df_marion_elements.metamer == phyto_id) & (df_marion_elements.organ == 'sheath')]['t'].unique(), df_marion_elements[(df_marion_elements.metamer == phyto_id) & (df_marion_elements.organ == 'sheath')].groupby('t')['N_tot'].sum(), color='orange', label='Control')
        axs[2][2].set_xlim(tmin, tmax)
        axs[2][2].set_ylabel('N content sh (g)')

        plt.tight_layout()
        pdf.savefig()  # saves the current figure into a pdf page
        plt.close()

def leaf_emergence(df_current_hz, df_marion_hz):
    dict_current, dict_marion = {'phyto_id': [], 't_emergence': []}, {'phyto_id': [], 't_emergence_marion': []}
    for phyto_id in df_current_hz['metamer'].unique():
        if not df_current_hz[(df_current_hz.metamer == phyto_id)]['leaf_is_emerged'].any():
            continue
        dict_current['phyto_id'].append(phyto_id)
        dict_current['t_emergence'].append(df_current_hz[(df_current_hz.metamer == phyto_id) & (df_current_hz.leaf_is_emerged == True)]['t'].iloc[0])
    for phyto_id in df_marion_hz['metamer'].unique():
        if phyto_id == 3 or not df_marion_hz[(df_marion_hz.metamer == phyto_id)]['leaf_is_emerged'].any():
            continue
        dict_marion['phyto_id'].append(phyto_id)
        dict_marion['t_emergence_marion'].append(df_marion_hz[(df_marion_hz.metamer == phyto_id) & (df_marion_hz.leaf_is_emerged == True)]['t'].iloc[0])

    df_current = pd.DataFrame.from_dict(dict_current)
    df_marion = pd.DataFrame.from_dict(dict_marion)
    df_merged = df_current.merge(df_marion, on='phyto_id', how='left')
    axis = df_merged.plot(kind='bar')
    axis.set_xlabel('N° de feuille')
    axis.set_ylabel('Temps emergence (hour)')

    ax2 = axis.twinx()
    ax2.set_yticks(axis.get_yticks())
    ax2.set_yticklabels(meteo_data.loc[axis.get_yticks()]['Date'].dt.strftime('%d/%m'))
    ax2.yaxis.set_ticks_position('left')  # set the position of the second x-axis to bottom
    ax2.yaxis.set_label_position('left')  # set the position of the second x-axis to bottom
    ax2.spines['left'].set_position(('outward', 50))

    fig = axis.get_figure()

    plt.tight_layout()
    pdf.savefig(fig)  # saves the current figure into a pdf page
    plt.close()

def plastochrone(df_current_axes, df_marion_axes):
    fig, axis = plt.subplots()

    axis.plot(df_current_axes['t'], df_current_axes['nb_leaves'], label='current')
    axis.plot(df_marion_axes['t'], df_marion_axes['nb_leaves'], label='Control')

    axis.set_xlabel('Time (day)')
    axis.set_ylabel('Number of leaves on MS')
    axis.legend()

    plt.tight_layout()
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()


if __name__ == '__main__':
    # get current directory
    # path = os.getcwd()
    
     ### Path Current

    #Simul CNWheat Plantfusion
    #path = r'C:\Users\agrumel\Code\Python_Ecophy\plantfusion_ag\outputs\cnwheat_default_temp\wheat'

    #simul CNWheat Marion depuis mon ordi
    #path = r'C:\Users\agrumel\Code\Python_Ecophy\WheatFspm\fspm-wheat\example\Vegetative_stages'

    #Simul CNWheat Soil3DS
    #path = r'C:\Users\agrumel\Code\Python_Ecophy\plantfusion_ag\outputs\cnwheat_soil3ds\bound\0.2m\wheat'
    #path = r'C:\Users\agrumel\Code\Python_Ecophy\plantfusion_ag\outputs\cnwheat_soil3ds\profile\0.2m\wheat'
    #path = r'C:\Users\agrumel\Code\Python_Ecophy\plantfusion_ag\outputs\cnwheat_soil3ds_debug\homogeneous\1.5m\wheat'
    #path = r'C:\Users\agrumel\Code\Python_Ecophy\plantfusion_ag\outputs\cnwheat_soil3ds_debug\profile\0.2m\wheat'

    #path = r'C:\Users\agrumel\Code\Python_Ecophy\plantfusion_ag\outputs\cnwheat_soil3ds\profile\0.2m_50N\wheat'
    path = r'C:\Users\agrumel\Code\Python_Ecophy\plantfusion_ag\outputs\cnwheat_soil3ds\profile\lowN\0.2m_50N\wheat'
    path = r'C:\Users\agrumel\Documents\Données\Sorties CNWheat\lowN_monoculm\0.2m_50N\wheat'
    path = r'C:\Users\agrumel\Documents\Données\Sorties CNWheat\lowN\0.2m_50N\wheat'


    #path = r'C:\Users\agrumel\Code\Python_Ecophy\plantfusion_ag\outputs\cnwheat_soil3ds\profile\54321_N\wheat'
    #path = r'C:\Users\agrumel\Code\Python_Ecophy\plantfusion_ag\outputs\cnwheat_soil3ds\profile\302010_N\wheat'

    #path = r'C:\Users\agrumel\Code\Python_Ecophy\plantfusion_ag\outputs\cnwheat_soil3ds\bound\0.2m_N\wheat'


    #Simul CNWheat tillers
    #path = r'C:\Users\agrumel\Code\Python_Ecophy\plantfusion_ag\outputs\cnwheat_default_tillers\0til\wheat'

    POSTPROCESSING = os.path.join(path, 'postprocessing')
    GRAPHS = os.path.join(path, 'graphs')

    ### Path Control
    
    #Version Soumission Marion
    #dirpath_control = r'C:\Users\agrumel\Documents\Données\Sorties CNWheat\Données Marion\Soumission_JXBot'

    #Simul CNWheat Plantfusion
    #dirpath_control = r'C:\Users\agrumel\Code\Python_Ecophy\plantfusion_ag\outputs\cnwheat_default\wheat'

    #simul CNWheat Marion depuis mon ordi
    #dirpath_control = r'C:\Users\agrumel\Code\Python_Ecophy\WheatFspm\fspm-wheat\example\Vegetative_stages'

    #Simul CNWheat OpenAlea
    #dirpath_control = r'C:\Users\agrumel\Documents\Données\Sorties CNWheat\Vegetative_stages - V2 Marion\outputs'


    OUTPUTS_MARION = os.path.join(dirpath_control, 'brut')
    POSTPROCESSING_MARION = os.path.join(dirpath_control, 'postprocessing')
    GRAPHS_MARION = os.path.join(dirpath_control, 'graphs')

    # Axes
    df_current_axes = pd.read_csv(os.path.join(POSTPROCESSING, 'axes_postprocessing.csv'))
    df_current_axes = df_current_axes[df_current_axes['axis'] == 'MS']
    df_marion_axes = pd.read_csv(os.path.join(POSTPROCESSING_MARION, 'axes_postprocessing.csv'))
    df_marion_axes = df_marion_axes[df_marion_axes['axis'] == 'MS']
    df_marion_axes['t'] = df_marion_axes['t'] + delta_t_simuls
    df_marion_axes_outputs = pd.read_csv(os.path.join(OUTPUTS_MARION, 'axes_outputs.csv'))
    df_marion_axes_outputs = df_marion_axes_outputs[df_marion_axes_outputs['axis'] == 'MS']
    df_marion_axes_outputs['t'] = df_marion_axes_outputs['t'] + delta_t_simuls

    # SAMs
    # df_marion_SAMS = pd.read_csv(os.path.join(dirpath_control, 'outputs', 'SAM_states.csv'))
    # df_marion_SAMS = df_marion_SAMS[df_marion_SAMS['axis'] == 'MS']
    # df_marion_SAMS['t'] = df_marion_SAMS['t'] + delta_t_simuls

    # Organs
    df_current_organs = pd.read_csv(os.path.join(POSTPROCESSING, 'organs_postprocessing.csv'))
    df_current_organs = df_current_organs[df_current_organs['axis'] == 'MS']
    df_marion_organs = pd.read_csv(os.path.join(POSTPROCESSING_MARION, 'organs_postprocessing.csv'))
    df_marion_organs = df_marion_organs[df_marion_organs['axis'] == 'MS']
    df_marion_organs['t'] = df_marion_organs['t'] + delta_t_simuls

    # Elements
    df_current_elements = pd.read_csv(os.path.join(POSTPROCESSING, 'elements_postprocessing.csv'))
    df_current_elements = df_current_elements[df_current_elements['axis'] == 'MS']
    df_marion_elements = pd.read_csv(os.path.join(POSTPROCESSING_MARION, 'elements_postprocessing.csv'))
    df_marion_elements = df_marion_elements[df_marion_elements['axis'] == 'MS']
    df_marion_elements['t'] = df_marion_elements['t'] + delta_t_simuls

    # HZ
    df_current_hz = pd.read_csv(os.path.join(POSTPROCESSING, 'hiddenzones_postprocessing.csv'))
    df_current_hz = df_current_hz[df_current_hz['axis'] == 'MS']
    df_marion_hz = pd.read_csv(os.path.join(POSTPROCESSING_MARION, 'hiddenzones_postprocessing.csv'))
    df_marion_hz = df_marion_hz[df_marion_hz['axis'] == 'MS']
    df_marion_hz['t'] = df_marion_hz['t'] + delta_t_simuls

    tmin = df_current_axes.t.min()
    tmax = df_current_axes.t.max()

    # plot graphs
    with PdfPages('Comparison_cnwheat.pdf') as pdf:
        # phloem
        phloem(df_current_organs, df_marion_organs, tmin, tmax)

        # Photosynthesis
        photosynthesis(df_current_axes, df_marion_axes)

        # roots
        roots(df_current_organs, df_marion_organs, tmin, tmax)

        # dry mass & shoot : root
        dry_mass(df_current_axes, df_marion_axes, df_current_organs, df_marion_organs, tmin, tmax)

        # N mass
        N_mass(df_current_axes, df_marion_axes, df_current_organs, df_marion_organs, tmin, tmax)

        # Surfaces
        surface(df_current_elements, df_marion_elements, tmin, tmax)
        include_images(GRAPHS, GRAPHS_MARION)

        # Leaf length & mstruct
        leaf_length_mstruct_area(df_current_hz, df_marion_hz, df_current_elements, df_marion_elements, tmin, tmax)

        # Leaf emergence date
        leaf_emergence(df_current_hz, df_marion_hz)

        # Plastochron
        plastochrone(df_current_axes, df_marion_axes_outputs)