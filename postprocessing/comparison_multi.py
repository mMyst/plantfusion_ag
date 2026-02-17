import os
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

delta_t_simuls = 0 #1509
meteo_data = pd.read_csv(os.path.join('inputs_fspmwheat', 'meteo_Ljutovac2002.csv'), index_col='t')
meteo_data['Date'] = pd.to_datetime(meteo_data['Date'], format='%d/%m/%Y')


def add_date_axis(ax, meteo_data):
    """Helper function to add the secondary Date x-axis."""
    ax2 = ax.twiny()
    ax2.set_xticks(ax.get_xticks())
    # Ensure ticks exist in meteo_data to prevent KeyError
    valid_ticks = [t for t in ax.get_xticks() if t in meteo_data.index]
    if valid_ticks:
        ax2.set_xticklabels(meteo_data.loc[valid_ticks]['Date'].dt.strftime('%d/%m'))
    ax2.xaxis.set_ticks_position('bottom')
    ax2.xaxis.set_label_position('bottom')
    ax2.spines['bottom'].set_position(('outward', 35))


def phloem(dict_current_organs, df_control_organs, tmin, tmax, pdf):
    fig, axs = plt.subplots(2, 2, figsize=(10, 8))

    # 1. Concentration sucrose
    axs[0, 0].plot(df_control_organs[df_control_organs.organ == 'phloem']['t'], df_control_organs[df_control_organs.organ == 'phloem']['Conc_Sucrose'], label='Control', color='black', linestyle='--')
    for name, df in dict_current_organs.items():
        axs[0, 0].plot(df[df.organ == 'phloem']['t'], df[df.organ == 'phloem']['Conc_Sucrose'], label=name)
    axs[0, 0].legend()
    axs[0, 0].set_xlim(tmin, tmax)
    axs[0, 0].set_ylabel('Concentration sucrose (µmol g-1)')

    # 2. Amount of sucrose
    axs[0, 1].plot(df_control_organs[df_control_organs.organ == 'phloem']['t'], df_control_organs[df_control_organs.organ == 'phloem']['sucrose'], label='Control', color='black', linestyle='--')
    for name, df in dict_current_organs.items():
        axs[0, 1].plot(df[df.organ == 'phloem']['t'], df[df.organ == 'phloem']['sucrose'], label=name)
    axs[0, 1].legend()
    axs[0, 1].set_xlim(tmin, tmax)
    axs[0, 1].set_ylabel('Amount of sucrose (µmol C)')

    # 3. Concentration amino acids
    axs[1, 0].plot(df_control_organs[df_control_organs.organ == 'phloem']['t'], df_control_organs[df_control_organs.organ == 'phloem']['Conc_Amino_Acids'], label='Control', color='black', linestyle='--')
    for name, df in dict_current_organs.items():
        axs[1, 0].plot(df[df.organ == 'phloem']['t'], df[df.organ == 'phloem']['Conc_Amino_Acids'], label=name)
    axs[1, 0].legend()
    axs[1, 0].set_xlim(tmin, tmax)
    axs[1, 0].set_ylabel('Concentration amino acids (µmol g-1)')
    add_date_axis(axs[1, 0], meteo_data)

    # 4. Amino acids amount
    axs[1, 1].plot(df_control_organs[df_control_organs.organ == 'phloem']['t'], df_control_organs[df_control_organs.organ == 'phloem']['amino_acids'], label='Control', color='black', linestyle='--')
    for name, df in dict_current_organs.items():
        axs[1, 1].plot(df[df.organ == 'phloem']['t'], df[df.organ == 'phloem']['amino_acids'], label=name)
    axs[1, 1].legend()
    axs[1, 1].set_xlim(tmin, tmax)
    axs[1, 1].set_ylabel('amino acids (µmol N)')
    add_date_axis(axs[1, 1], meteo_data)

    plt.tight_layout()
    pdf.savefig()
    plt.close()


def photosynthesis(dict_current_axes, df_control_axes, pdf):
    fig, axis = plt.subplots()
    
    df_control_axes['day'] = df_control_axes['t'] // 24 + 1
    axis.plot(df_control_axes['day'].unique(), df_control_axes.groupby('day')['Total_Photosynthesis'].sum(), label='Control', color='black', linestyle='--')
    
    for name, df in dict_current_axes.items():
        df['day'] = df['t'] // 24 + 1
        axis.plot(df['day'].unique(), df.groupby('day')['Total_Photosynthesis'].sum(), label=name)

    axis.set_xlabel('Time (day)')
    axis.set_ylabel('Total Photosynthesis µmol C')
    axis.legend()

    plt.tight_layout()
    pdf.savefig()
    plt.close()


def roots(dict_current_organs, df_control_organs, tmin, tmax, pdf):
    fig, axs = plt.subplots(2, 2, figsize=(10, 8))

    # 1. Concentration sucrose
    axs[0, 0].plot(df_control_organs[df_control_organs.organ == 'roots']['t'], df_control_organs[df_control_organs.organ == 'roots']['Conc_Sucrose'], label='Control', color='black', linestyle='--')
    for name, df in dict_current_organs.items():
        axs[0, 0].plot(df[df.organ == 'roots']['t'], df[df.organ == 'roots']['Conc_Sucrose'], label=name)
    axs[0, 0].legend()
    axs[0, 0].set_xlim(tmin, tmax)
    axs[0, 0].set_ylabel('Concentration sucrose (µmol g-1)')

    # 2. Amount of sucrose
    axs[0, 1].plot(df_control_organs[df_control_organs.organ == 'roots']['t'], df_control_organs[df_control_organs.organ == 'roots']['sucrose'], label='Control', color='black', linestyle='--')
    for name, df in dict_current_organs.items():
        axs[0, 1].plot(df[df.organ == 'roots']['t'], df[df.organ == 'roots']['sucrose'], label=name)
    axs[0, 1].legend()
    axs[0, 1].set_xlim(tmin, tmax)
    axs[0, 1].set_ylim(0, 1000)
    axs[0, 1].set_ylabel('Amount of sucrose (µmol C)')

    # 3. Concentration nitrates
    axs[1, 0].plot(df_control_organs[df_control_organs.organ == 'roots']['t'], df_control_organs[df_control_organs.organ == 'roots']['Conc_Nitrates'], label='Control', color='black', linestyle='--')
    for name, df in dict_current_organs.items():
        axs[1, 0].plot(df[df.organ == 'roots']['t'], df[df.organ == 'roots']['Conc_Nitrates'], label=name)
    axs[1, 0].legend()
    axs[1, 0].set_xlim(tmin, tmax)
    axs[1, 0].set_ylabel('Concentration nitrates (µmol g-1)')
    add_date_axis(axs[1, 0], meteo_data)

    # 4. Cytokinins
    axs[1, 1].plot(df_control_organs[df_control_organs.organ == 'roots']['t'], df_control_organs[df_control_organs.organ == 'roots']['Conc_cytokinins'], label='Control', color='black', linestyle='--')
    for name, df in dict_current_organs.items():
        axs[1, 1].plot(df[df.organ == 'roots']['t'], df[df.organ == 'roots']['Conc_cytokinins'], label=name)
    axs[1, 1].legend()
    axs[1, 1].set_xlim(tmin, tmax)
    axs[1, 1].set_ylabel('Conc_cytokinins (AU g-1)')
    add_date_axis(axs[1, 1], meteo_data)

    plt.tight_layout()
    pdf.savefig()
    plt.close()


def dry_mass(dict_current_axes, df_control_axes, dict_current_organs, df_control_organs, tmin, tmax, pdf):
    fig, axs = plt.subplots(2, 2, sharex=True, figsize=(10, 8))

    # Dry mass shoot
    axs[0, 0].plot(df_control_axes['t'], df_control_axes['sum_dry_mass_shoot'], label='Control', color='black', linestyle='--')
    for name, df in dict_current_axes.items():
        axs[0, 0].plot(df['t'], df['sum_dry_mass_shoot'], label=name)
    axs[0, 0].legend()
    axs[0, 0].set_xlim(tmin, tmax)
    axs[0, 0].set_ylabel('Dry mass shoot (g)')

    # Dry mass roots
    axs[0, 1].plot(df_control_axes['t'], df_control_axes['sum_dry_mass_roots'], label='Control', color='black', linestyle='--')
    for name, df in dict_current_axes.items():
        axs[0, 1].plot(df['t'], df['sum_dry_mass_roots'], label=name)
    axs[0, 1].legend()
    axs[0, 1].set_xlim(tmin, tmax)
    axs[0, 1].set_ylabel('Dry mass roots (g)')

    # mstruct shoot
    axs[1, 0].plot(df_control_axes['t'], df_control_axes['mstruct_shoot'], label='Control', color='black', linestyle='--')
    for name, df in dict_current_axes.items():
        axs[1, 0].plot(df['t'], df['mstruct_shoot'], label=name)
    axs[1, 0].legend()
    axs[1, 0].set_xlim(tmin, tmax)
    axs[1, 0].set_ylabel('mstruct shoot (g)')
    add_date_axis(axs[1, 0], meteo_data)

    # mstruct roots
    axs[1, 1].plot(df_control_organs[df_control_organs['organ'] == 'roots']['t'], df_control_organs[df_control_organs['organ'] == 'roots']['mstruct'], label='Control', color='black', linestyle='--')
    for name, df in dict_current_organs.items():
        axs[1, 1].plot(df[df['organ'] == 'roots']['t'], df[df['organ'] == 'roots']['mstruct'], label=name)
    axs[1, 1].legend()
    axs[1, 1].set_xlim(tmin, tmax)
    axs[1, 1].set_ylabel('mstruct roots (g)')
    add_date_axis(axs[1, 1], meteo_data)

    plt.tight_layout()
    pdf.savefig()
    plt.close()

    # shoot : root
    fig, axis = plt.subplots()
    axis.plot(df_control_axes['t'], df_control_axes['shoot_roots_ratio'], label='Control', color='black', linestyle='--')
    for name, df in dict_current_axes.items():
        axis.plot(df['t'], df['shoot_roots_ratio'], label=name)
    axis.legend()
    axis.set_xlim(tmin, tmax)
    axis.set_ylim(0, 2)
    axis.set_ylabel('shoot : root ratio')
    add_date_axis(axis, meteo_data)
    
    plt.tight_layout()
    pdf.savefig()
    plt.close()


def N_mass(dict_current_axes, df_control_axes, dict_current_organs, df_control_organs, tmin, tmax, pdf):
    fig, axs = plt.subplots(2, 2, sharex=True, figsize=(10, 8))

    # % N axis
    axs[0, 0].plot(df_control_axes['t'], df_control_axes['N_content'], label='Control', color='black', linestyle='--')
    for name, df in dict_current_axes.items():
        axs[0, 0].plot(df['t'], df['N_content'], label=name)
    axs[0, 0].legend()
    axs[0, 0].set_xlim(tmin, tmax)
    axs[0, 0].set_ylim(0, 10)
    axs[0, 0].set_ylabel('N content axis (% DM)')

    # N shoot
    axs[0, 1].plot(df_control_axes['t'], df_control_axes['N_content_shoot'], label='Control', color='black', linestyle='--')
    for name, df in dict_current_axes.items():
        axs[0, 1].plot(df['t'], df['N_content_shoot'], label=name)
    axs[0, 1].legend()
    axs[0, 1].set_xlim(tmin, tmax)
    axs[0, 1].set_ylabel('N content shoot (% DM)')

    # N axis
    axs[1, 0].plot(df_control_axes['t'], df_control_axes['sum_N_g'], label='Control', color='black', linestyle='--')
    for name, df in dict_current_axes.items():
        axs[1, 0].plot(df['t'], df['sum_N_g'], label=name)
    axs[1, 0].legend()
    axs[1, 0].set_xlim(tmin, tmax)
    axs[1, 0].set_ylabel('N content axis (g)')
    add_date_axis(axs[1, 0], meteo_data)

    # N uptake
    axs[1, 1].plot(df_control_organs[df_control_organs['organ'] == 'roots']['t'], df_control_organs[df_control_organs['organ'] == 'roots']['Uptake_Nitrates'], label='Control', color='black', linestyle='--')
    for name, df in dict_current_organs.items():
        axs[1, 1].plot(df[df['organ'] == 'roots']['t'], df[df['organ'] == 'roots']['Uptake_Nitrates'], label=name)
    axs[1, 1].legend()
    axs[1, 1].set_xlim(tmin, tmax)
    axs[1, 1].set_ylabel('Nitrate uptake (µmol)')
    add_date_axis(axs[1, 1], meteo_data)

    plt.tight_layout()
    pdf.savefig()
    plt.close()


def surface(dict_current_elements, df_control_elements, tmin, tmax, pdf):
    fig, axs = plt.subplots(2, 2, figsize=(10, 8))

    # Total green area
    axs[0, 0].plot(df_control_elements['t'].unique(), df_control_elements.groupby('t')['green_area'].sum(), label='Control', color='black', linestyle='--')
    for name, df in dict_current_elements.items():
        axs[0, 0].plot(df['t'].unique(), df.groupby('t')['green_area'].sum(), label=name)
    axs[0, 0].legend()
    axs[0, 0].set_xlim(tmin, tmax)
    axs[0, 0].set_ylabel('Total green area (m²)')

    # Blade green area
    df_control_elements_blade = df_control_elements[df_control_elements.organ == 'blade']
    axs[0, 1].plot(df_control_elements_blade['t'].unique(), df_control_elements_blade.groupby('t')['green_area'].sum(), label='Control', color='black', linestyle='--')
    for name, df in dict_current_elements.items():
        df_blade = df[df.organ == 'blade']
        axs[0, 1].plot(df_blade['t'].unique(), df_blade.groupby('t')['green_area'].sum(), label=name)
    axs[0, 1].legend()
    axs[0, 1].set_xlim(tmin, tmax)
    axs[0, 1].set_ylabel('Blade green area (m²)')

    # Sheath green area
    df_control_elements_sheath = df_control_elements[df_control_elements.organ == 'sheath']
    axs[1, 0].plot(df_control_elements_sheath['t'].unique(), df_control_elements_sheath.groupby('t')['green_area'].sum(), label='Control', color='black', linestyle='--')
    for name, df in dict_current_elements.items():
        df_sheath = df[df.organ == 'sheath']
        axs[1, 0].plot(df_sheath['t'].unique(), df_sheath.groupby('t')['green_area'].sum(), label=name)
    axs[1, 0].legend()
    axs[1, 0].set_xlim(tmin, tmax)
    axs[1, 0].set_ylabel('Sheath green area (m²)')
    add_date_axis(axs[1, 0], meteo_data)

    # Internode green area
    df_control_elements_internode = df_control_elements[df_control_elements.organ == 'internode']
    axs[1, 1].plot(df_control_elements_internode['t'].unique(), df_control_elements_internode.groupby('t')['green_area'].sum(), label='Control', color='black', linestyle='--')
    for name, df in dict_current_elements.items():
        df_internode = df[df.organ == 'internode']
        axs[1, 1].plot(df_internode['t'].unique(), df_internode.groupby('t')['green_area'].sum(), label=name)
    axs[1, 1].legend()
    axs[1, 1].set_xlim(tmin, tmax)
    axs[1, 1].set_ylabel('Internode green area (m²)')
    add_date_axis(axs[1, 1], meteo_data)

    plt.tight_layout()
    pdf.savefig()
    plt.close()


def height(dict_current_elements, df_control_elements, tmin, tmax, pdf):
    fig, axis = plt.subplots()

    # Plant height
    axis.plot(df_control_elements['t'].unique(), df_control_elements.groupby('t')['height'].max(), label='Control', color='black', linestyle='--')
    for name, df in dict_current_elements.items():
        axis.plot(df['t'].unique(), df.groupby('t')['height'].max(), label=name)
    axis.legend()
    axis.set_xlim(tmin, tmax)
    axis.set_ylabel('Plant Height (m)')
    add_date_axis(axis, meteo_data)

    plt.tight_layout()
    pdf.savefig()
    plt.close()


def include_images(dict_graphs_paths, graphs_Marion_path, pdf):
    graphs_to_include = ['leaf_L_hz.PNG', 'Leaf_Lmax.PNG', 'RER_comparison.PNG', 'phyllochron.PNG', 'lamina_Wmax.PNG', 'SSLW.PNG']
    for graph in graphs_to_include:
        # Loop through simulations to display their graphs iteratively
        for sim_name, graphs_path in dict_graphs_paths.items():
            img_path = os.path.join(graphs_path, graph)
            if os.path.exists(img_path):
                im = plt.imread(img_path)
                fig = plt.figure(figsize=(13, 10))
                fig.figimage(im)
                fig.suptitle(f"{sim_name} - {graph}")
                pdf.savefig()
                plt.close()


def leaf_length_mstruct_area(dict_current_hz, df_control_hz, dict_current_elements, df_control_elements, tmin, tmax, pdf):
    # Loop through phytomer using Control as the base to define metamers
    for phyto_id in df_control_hz['metamer'].unique():
        fig, axs = plt.subplots(3, 3, figsize=(12, 10))

        # 1. Length leaf
        if phyto_id in (1, 2):
            axs[0][0].plot(df_control_elements[(df_control_elements.metamer == phyto_id)]['t'].unique(),
                            df_control_elements[(df_control_elements.metamer == phyto_id)].groupby('t')['length'].sum(), 
                            label='Control', color='black', linestyle='--')
        else:
            axs[0][0].plot(df_control_hz[(df_control_hz.metamer == phyto_id)]['t'].unique(),
                            df_control_hz[(df_control_hz.metamer == phyto_id)].groupby('t')['leaf_L'].sum(), 
                            label='Control', color='black', linestyle='--')
        for name, df in dict_current_hz.items():
            axs[0][0].plot(df[(df.metamer == phyto_id)]['t'].unique(), 
                           df[(df.metamer == phyto_id)].groupby('t')['leaf_L'].sum(), 
                           label=name)
        axs[0][0].set_xlim(tmin, tmax)
        axs[0][0].set_ylabel('Leaf_L  hz' + ' (m)')
        axs[0][0].set_title('Leaf' + '_' + str(phyto_id))
        axs[0][0].set_xticks([])

        # 2. Length sheath
        axs[0][1].plot(df_control_elements[(df_control_elements.metamer == phyto_id) & (df_control_elements.organ == 'sheath')]['t'].unique(), df_control_elements[(df_control_elements.metamer == phyto_id) & (df_control_elements.organ == 'sheath')].groupby('t')['length'].sum(), label='Control', color='black', linestyle='--')
        for name, df in dict_current_elements.items():
            axs[0][1].plot(df[(df.metamer == phyto_id) & (df.organ == 'sheath')]['t'].unique(),
                            df[(df.metamer == phyto_id) & (df.organ == 'sheath')].groupby('t')['length'].sum(),
                              label=name)
        axs[0][1].set_xlim(tmin, tmax)
        axs[0][1].set_ylabel('Sheath_L elt' + ' (m)')
        axs[0][1].set_xticks([])

        # 3. Internode length
        Marion_internode_L = df_control_hz[df_control_hz.metamer == phyto_id].groupby('t')['internode_L'].sum()
        axs[0][2].plot(df_control_hz[(df_control_hz.metamer == phyto_id)]['t'].unique(), 
                       Marion_internode_L,
                        label='Control', color='black', linestyle='--')
        for name, df in dict_current_hz.items():
            axs[0][2].plot(df[(df.metamer == phyto_id)]['t'].unique(),
                            df[df.metamer == phyto_id].groupby('t')['internode_L'].sum(),
                            label=name)
        axs[0][2].set_xlim(tmin, tmax)
        axs[0][2].set_ylabel('Internode_L hz (m)')
        axs[0][2].legend(loc='upper center', bbox_to_anchor=(0.25, 1.5), ncol=2, fontsize="8")

        # 4. Mstruct hz
        if phyto_id in (1, 2):
            axs[1][0].plot(df_control_elements[(df_control_elements.metamer == phyto_id)]['t'].unique(),
                          df_control_elements[(df_control_elements.metamer == phyto_id) & (df_control_elements.element == 'LeafElement1')]
                          .groupby(['t'])['mstruct'].sum(),
                          label='Control', color='black', linestyle='--')
            for name, df in dict_current_elements.items():
                axs[1][0].plot(df[(df.metamer == phyto_id) & (df.element == 'LeafElement1')]['t'].unique(), 
                               df[(df.metamer == phyto_id) & (df.element == 'LeafElement1')]
                               .groupby(['t'])['mstruct'].sum(), 
                               label=name)
        else:
            axs[1][0].plot(df_control_hz[(df_control_hz.metamer == phyto_id)]['t'].unique(),
                            df_control_hz[(df_control_hz.metamer == phyto_id)]
                            .groupby('t')['mstruct'].sum(), 
                            label='Control', color='black', linestyle='--')
            for name, df in dict_current_hz.items():
                axs[1][0].plot(df[(df.metamer == phyto_id)]['t'].unique(),
                                df[(df.metamer == phyto_id)]
                                .groupby('t')['mstruct'].sum(), 
                                label=name)
        axs[1][0].set_xlim(tmin, tmax)
        axs[1][0].set_ylabel('mstruct hz' + ' (g)')
        axs[1][0].set_xticks([])

        # 5. mstruct lamina
        if phyto_id != 0:
            axs[1][1].plot(df_control_elements[(df_control_elements.metamer == phyto_id) & (df_control_elements.organ == 'blade')]['t'].unique(),
                            df_control_elements[(df_control_elements.metamer == phyto_id) & (df_control_elements.organ == 'blade')]
                            .groupby(['t'])['mstruct'].sum(), 
                            label='Control', color='black', linestyle='--')
            for name, df in dict_current_elements.items():
                axs[1][1].plot(df[(df.metamer == phyto_id) & (df.organ == 'blade')]['t'].unique(),
                                df[(df.metamer == phyto_id) & (df.organ == 'blade')]
                                .groupby('t')['mstruct'].sum(), 
                                label=name)
            axs[1][1].set_xlim(tmin, tmax)
            axs[1][1].set_ylabel('Lam mstruct elt' + ' (g)')
            axs[1][1].set_xticks([])

        # 6. N content lamina
        if phyto_id != 0:
            Marion_lamina_L = df_control_elements[(df_control_elements.metamer == phyto_id) & (df_control_elements.element == 'LeafElement1')].groupby(['t'])['N_tot'].sum()
            axs[1][2].plot(df_control_elements[(df_control_elements.metamer == phyto_id) & (df_control_elements.element == 'LeafElement1')]['t'].unique(),
                            Marion_lamina_L, 
                            label='Control', color='black', linestyle='--')
            for name, df in dict_current_elements.items():
                axs[1][2].plot(df[(df.metamer == phyto_id) & (df.element == 'LeafElement1')]['t'].unique(),
                                df[(df.metamer == phyto_id) & (df.element == 'LeafElement1')]
                                .groupby(['t'])['N_tot'].sum(),
                                  label=name)
            axs[1][2].set_xlim(tmin, tmax)
            axs[1][2].set_ylabel('N content lam (g)')

        # 7. mstruct sheath
        axs[2][0].plot(df_control_elements[(df_control_elements.metamer == phyto_id) & (df_control_elements.organ == 'sheath')]['t'].unique(),
                        df_control_elements[(df_control_elements.metamer == phyto_id) & (df_control_elements.organ == 'sheath')]
                        .groupby('t')['mstruct'].sum(), 
                        label='Control', color='black', linestyle='--')
        for name, df in dict_current_elements.items():
            axs[2][0].plot(df[(df.metamer == phyto_id) & (df.organ == 'sheath')]['t'].unique(),
                            df[(df.metamer == phyto_id) & (df.organ == 'sheath')]
                            .groupby('t')['mstruct'].sum(), 
                            label=name)
        axs[2][0].set_xlim(tmin, tmax)
        axs[2][0].set_ylabel('Sheath mstruct' + ' (g)')

        # 8. Green area lamina
        Marion_blade_green_area = df_control_elements[(df_control_elements.metamer == phyto_id) & (df_control_elements.element == 'LeafElement1')]
        axs[2][1].plot(Marion_blade_green_area['t'].unique(),
                        Marion_blade_green_area.groupby('t')['green_area'].sum(),
                        label='Control', color='black', linestyle='--')
        for name, df in dict_current_elements.items():
            df_blade_ga = df[(df.metamer == phyto_id) & (df.element == 'LeafElement1')]
            axs[2][1].plot(df_blade_ga['t'].unique(),
                            df_blade_ga.groupby('t')['green_area'].sum(),
                            label=name)
        axs[2][1].set_xlim(tmin, tmax)
        axs[2][1].set_ylabel('Lamina GA' + ' (m2)')

        # 9. N content sheath
        if phyto_id != 0:
            axs[2][2].plot(df_control_elements[(df_control_elements.metamer == phyto_id) & (df_control_elements.organ == 'sheath')]['t'].unique(), df_control_elements[(df_control_elements.metamer == phyto_id) & (df_control_elements.organ == 'sheath')].groupby('t')['N_tot'].sum(), label='Control', color='black', linestyle='--')
            for name, df in dict_current_elements.items():
                axs[2][2].plot(df[(df.metamer == phyto_id) & (df.organ == 'sheath')]['t'].unique(), df[(df.metamer == phyto_id) & (df.organ == 'sheath')].groupby('t')['N_tot'].sum(), label=name)
            axs[2][2].set_xlim(tmin, tmax)
            axs[2][2].set_ylabel('N content sh (g)')

        plt.tight_layout()
        pdf.savefig()
        plt.close()


def leaf_emergence(dict_current_hz, df_control_hz, pdf):
    dict_all_emergence = {}

    # Gather Control Emergence
    t_emergence_marion = {}
    for phyto_id in df_control_hz['metamer'].unique():
        if phyto_id == 3 or not df_control_hz[(df_control_hz.metamer == phyto_id)]['leaf_is_emerged'].any():
            continue
        t_emergence_marion[phyto_id] = df_control_hz[(df_control_hz.metamer == phyto_id) & (df_control_hz.leaf_is_emerged == True)]['t'].iloc[0]
    dict_all_emergence['Control'] = t_emergence_marion

    # Gather Current Simulations Emergence
    for name, df in dict_current_hz.items():
        t_emergence = {}
        for phyto_id in df['metamer'].unique():
            if not df[(df.metamer == phyto_id)]['leaf_is_emerged'].any():
                continue
            t_emergence[phyto_id] = df[(df.metamer == phyto_id) & (df.leaf_is_emerged == True)]['t'].iloc[0]
        dict_all_emergence[name] = t_emergence

    # Plot as a grouped bar chart
    df_merged = pd.DataFrame(dict_all_emergence)
    axis = df_merged.plot(kind='bar', figsize=(10, 6))
    axis.set_xlabel('N° de feuille')
    axis.set_ylabel('Temps emergence (hour)')

    ax2 = axis.twinx()
    ax2.set_yticks(axis.get_yticks())
    
    valid_ticks = [t for t in axis.get_yticks() if t in meteo_data.index]
    if valid_ticks:
        ax2.set_yticklabels(meteo_data.loc[valid_ticks]['Date'].dt.strftime('%d/%m'))
        
    ax2.yaxis.set_ticks_position('left')
    ax2.yaxis.set_label_position('left')
    ax2.spines['left'].set_position(('outward', 50))

    fig = axis.get_figure()
    plt.tight_layout()
    pdf.savefig(fig)
    plt.close()


def plastochrone(dict_current_axes, df_control_axes, pdf):
    fig, axis = plt.subplots()

    axis.plot(df_control_axes['t'], df_control_axes['nb_leaves'], label='Control', color='black', linestyle='--')
    for name, df in dict_current_axes.items():
        axis.plot(df['t'], df['nb_leaves'], label=name)

    axis.set_xlabel('Time (day)')
    axis.set_ylabel('Number of leaves on MS')
    axis.legend()

    plt.tight_layout()
    pdf.savefig()
    plt.close()


if __name__ == '__main__':
    # 1. Define your dictionary of simulation paths

    simulations_paths = {
        # 'Sim_Default': r'C:\Users\agrumel\Code\Python_Ecophy\plantfusion_ag\outputs\cnwheat_default_temp\wheat',
        
        'Sim_Tillers_0': r'C:\Users\agrumel\Code\Python_Ecophy\plantfusion_ag\outputs\cnwheat_default_tillers\0til\wheat',
        # 'Sim_Tillers_1': r'C:\Users\agrumel\Code\Python_Ecophy\plantfusion_ag\outputs\cnwheat_default_tillers\1til\wheat',
        # 'Sim_Tillers_2': r'C:\Users\agrumel\Code\Python_Ecophy\plantfusion_ag\outputs\cnwheat_default_tillers\2til\wheat',
        # 'Sim_Tillers_3': r'C:\Users\agrumel\Code\Python_Ecophy\plantfusion_ag\outputs\cnwheat_default_tillers\3til\wheat'
        # Add as many as needed here...

        #'Sim_Soil3DS_Bound_0.2m': r'C:\Users\agrumel\Code\Python_Ecophy\plantfusion_ag\outputs\cnwheat_soil3ds\bound\0.2m\wheat',
        # 'Sim_Soil3DS_Profile_0.2m': r'C:\Users\agrumel\Code\Python_Ecophy\plantfusion_ag\outputs\cnwheat_soil3ds\profile\0.2m\wheat',
        # 'Sim_Soil3DS_Debug_Homogeneous_1.5m': r'C:\Users\agrumel\Code\Python_Ecophy\plantfusion_ag\outputs\cnwheat_soil3ds_debug\homogeneous\1.5m\wheat',
        # 'Sim_Soil3DS_Debug_Profile_0.2m': r'C:\Users\agrumel\Code\Python_Ecophy\plantfusion_ag\outputs\cnwheat_soil3ds_debug\profile\0.2m\wheat'
 
        # 'Sim_Soil3DS_Profile_0.2m_50N': r'C:\Users\agrumel\Code\Python_Ecophy\plantfusion_ag\outputs\cnwheat_soil3ds\profile\0.2m_50N\wheat',
        # 'Sim_Soil3DS_Profile_lowN_0.2m_50N': r'C:\Users\agrumel\Code\Python_Ecophy\plantfusion_ag\outputs\cnwheat_soil3ds\profile\lowN\0.2m_50N\wheat',
        # 'Sim_lowN_monoculm_0.2m_50N': r'C:\Users\agrumel\Documents\Données\Sorties CNWheat\lowN_monoculm\0.2m_50N\wheat',
        # 'Sim_lowN_0.2m_50N': r'C:\Users\agrumel\Documents\Données\Sorties CNWheat\lowN\0.2m_50N\wheat'




    }

    ### Path Control
    
    #Version Soumission Marion
    #dirpath_control = r'C:\Users\agrumel\Documents\Données\Sorties CNWheat\Données Marion\Soumission_JXBot'

    #Simul CNWheat Plantfusion
    dirpath_control = r'C:\Users\agrumel\Code\Python_Ecophy\plantfusion_ag\outputs\cnwheat_default\wheat'

    #simul CNWheat Marion depuis mon ordi
    #dirpath_control = r'C:\Users\agrumel\Code\Python_Ecophy\WheatFspm\fspm-wheat\example\Vegetative_stages'

    #Simul CNWheat OpenAlea
    #dirpath_control = r'C:\Users\agrumel\Documents\Données\Sorties CNWheat\Vegetative_stages - V2 Marion\outputs'

    

    OUTPUTS_CONTROL = os.path.join(dirpath_control, 'brut')
    POSTPROCESSING_CONTROL = os.path.join(dirpath_control, 'postprocessing')
    GRAPHS_CONTROL = os.path.join(dirpath_control, 'graphs')

    # Load Control Data
    df_control_axes = pd.read_csv(os.path.join(POSTPROCESSING_CONTROL, 'axes_postprocessing.csv'))
    #df_control_axes = df_control_afxes[df_control_axes['axis'] == 'MS']
    df_control_axes['t'] = df_control_axes['t'] + delta_t_simuls
    
    df_control_axes_outputs = pd.read_csv(os.path.join(OUTPUTS_CONTROL, 'axes_outputs.csv'))
    #df_control_axes_outputs = df_control_axes_outputs[df_control_axes_outputs['axis'] == 'MS']
    df_control_axes_outputs['t'] = df_control_axes_outputs['t'] + delta_t_simuls

    df_control_organs = pd.read_csv(os.path.join(POSTPROCESSING_CONTROL, 'organs_postprocessing.csv'))
    #df_control_organs = df_control_organs[df_control_organs['axis'] == 'MS']
    df_control_organs['t'] = df_control_organs['t'] + delta_t_simuls

    df_control_elements = pd.read_csv(os.path.join(POSTPROCESSING_CONTROL, 'elements_postprocessing.csv'))
    #df_control_elements = df_control_elements[df_control_elements['axis'] == 'MS']
    df_control_elements['t'] = df_control_elements['t'] + delta_t_simuls

    df_control_hz = pd.read_csv(os.path.join(POSTPROCESSING_CONTROL, 'hiddenzones_postprocessing.csv'))
    #df_control_hz = df_control_hz[df_control_hz['axis'] == 'MS']
    df_control_hz['t'] = df_control_hz['t'] + delta_t_simuls

    # 2. Load Current Simulations Data
    dict_current_axes = {}
    dict_current_organs = {}
    dict_current_elements = {}
    dict_current_hz = {}
    dict_graphs_paths = {}
    dict_current_axes_outputs = {}

    for sim_name, path in simulations_paths.items():
        POSTPROCESSING = os.path.join(path, 'postprocessing')
        dict_graphs_paths[sim_name] = os.path.join(path, 'graphs')

        # Load and filter data for each simulation
        df_ax = pd.read_csv(os.path.join(POSTPROCESSING, 'axes_postprocessing.csv'))
        df_org = pd.read_csv(os.path.join(POSTPROCESSING, 'organs_postprocessing.csv'))
        df_elt = pd.read_csv(os.path.join(POSTPROCESSING, 'elements_postprocessing.csv'))
        df_hz = pd.read_csv(os.path.join(POSTPROCESSING, 'hiddenzones_postprocessing.csv'))

        # Filter for main stem (MS) only
        # df_ax= df_ax[df_ax['axis'] == 'MS']
        # df_org = df_org[df_org['axis'] == 'MS']
        # df_elt = df_elt[df_elt['axis'] == 'MS']
        # df_hz = df_hz[df_hz['axis'] == 'MS']

        # Update dictionaries with filtered data
        dict_current_axes[sim_name] = df_ax
        dict_current_organs[sim_name] = df_org
        dict_current_elements[sim_name] = df_elt
        dict_current_hz[sim_name] = df_hz

        OUTPUTS = os.path.join(path, 'brut')
        df_ax_outputs = pd.read_csv(os.path.join(OUTPUTS, 'axes_outputs.csv'))
        #df_ax_outputs = df_ax_outputs[df_ax_outputs['axis'] == 'MS']
        dict_current_axes_outputs[sim_name] =  df_ax_outputs
        
        
    # Create filtered dictionaries and dataframes for MS only
    dict_current_hz_MS = {key: df[df['axis'] == "MS"] for key, df in dict_current_hz.items()}
    df_control_hz_MS = df_control_hz[df_control_hz['axis'] == 'MS']
    dict_current_elements_MS = {key: df[df['axis'] == "MS"] for key, df in dict_current_elements.items()}
    df_control_elements_MS = df_control_elements[df_control_elements['axis'] == 'MS']
    dict_current_axes_outputs_MS = {key: df[df['axis'] == "MS"] for key, df in dict_current_axes_outputs.items()}
    df_control_axes_outputs_MS = df_control_axes_outputs[df_control_axes_outputs['axis'] == 'MS']


    # Set tmin/tmax based on the first simulation available
    first_sim_key = list(dict_current_axes.keys())[0]
    tmin = dict_current_axes[first_sim_key].t.min()
    tmax = dict_current_axes[first_sim_key].t.max()

    # 3. Plot Graphs into PDF
    with PdfPages('Comparison_cnwheat.pdf') as pdf:
        phloem(dict_current_organs, df_control_organs, tmin, tmax, pdf)
        photosynthesis(dict_current_axes, df_control_axes, pdf)
        roots(dict_current_organs, df_control_organs, tmin, tmax, pdf)
        dry_mass(dict_current_axes, df_control_axes, dict_current_organs, df_control_organs, tmin, tmax, pdf)
        N_mass(dict_current_axes, df_control_axes, dict_current_organs, df_control_organs, tmin, tmax, pdf)
        surface(dict_current_elements, df_control_elements, tmin, tmax, pdf)
        height(dict_current_elements, df_control_elements, tmin, tmax, pdf)
        include_images(dict_graphs_paths, GRAPHS_CONTROL, pdf)
        leaf_length_mstruct_area(dict_current_hz_MS, df_control_hz_MS, dict_current_elements_MS, df_control_elements_MS, tmin, tmax, pdf)
        leaf_emergence(dict_current_hz_MS, df_control_hz_MS, pdf)
        plastochrone(dict_current_axes_outputs_MS, df_control_axes_outputs_MS, pdf)
