from plantfusion.new_wheat_wrapper import Wheat_wrapper
from plantfusion.l_egume_wrapper import L_egume_wrapper
from plantfusion.light_wrapper import Light_wrapper
from plantfusion.soil_wrapper import Soil_wrapper
from plantfusion.indexer import Indexer
from plantfusion.planter import Planter, grid_pattern_per_cell
from elongwheat import tillering as elongwheat_tillering

import numpy as np
import time
import datetime
import os
import math
import concurrent.futures

def simulation(
    in_folder_legume, in_folder_wheat, out_folder,
    start_wheat, simulation_length, id_usm, 
    run_postprocessing=False, run_graphs=False, writegeo=False, geostep=1,
    wheatroot_type='profile', min_depth=None
):
    try:
        # Create target Directory
        os.mkdir(os.path.normpath(out_folder))
        print("Directory ", os.path.normpath(out_folder), " Created ")
    except FileExistsError:
        print("Directory ", os.path.normpath(out_folder), " already exists")

    ######################
    ### INITIALIZATION ###
    ######################

    base_wheat_name = "wheat"
    legume_name = "legume"

    tillers_replications = {"T1": 0.5, "T2": 0.5, "T3": 0.5, "T4": 0.5}
    sky = "turtle46"
    RERmax_vegetative_stages_example = {
        "elongwheat": {
            "RERmax": {5: 3.35e-06, 6: 2.1e-06, 7: 2.0e-06, 8: 1.83e-06, 9: 1.8e-06, 10: 1.65e-06, 11: 1.56e-06}
        }
    }
    senescwheat_timestep = 1
    light_timestep = 4
  
    ### GENERATING POSITIONS
    col_pattern = (legume_name, "inter_row", base_wheat_name, base_wheat_name, "inter_row", legume_name)
  
    n_rows = 6
    n_cols = 6
    cell_size = 0.05 # Taille de la cellule en mètres
    
    offset = {base_wheat_name: 0.15*2, legume_name: 0}
    noise = {base_wheat_name:0.001, legume_name:0.01}

    # Utilisation de grid_pattern_per_cell pour individualiser les blés
    domain, full_grid_per_cell = grid_pattern_per_cell(
        n_cols=n_cols, 
        n_rows=n_rows, 
        cell_size=cell_size, 
        col_pattern=col_pattern, 
        noise=noise
    )

    # Filtrage : On récupère les blés individuels ("wheat_0", "wheat_1", etc.)
    wheat_names = [k for k in full_grid_per_cell.keys() if k.startswith(f"{base_wheat_name}_")]
    
    # Regroupement : On rassemble toutes les positions de luzerne sous une seule clé
    legume_positions = []
    for k in full_grid_per_cell.keys():
        if k.startswith(f"{legume_name}_"):
            legume_positions.extend(full_grid_per_cell[k])
            
    positions_grid = {legume_name: legume_positions}
    for w_name in wheat_names:
        positions_grid[w_name] = full_grid_per_cell[w_name]

    # Indexer: on passe le legume global, plus chaque instance de blé
    indexer = Indexer(global_order=[legume_name] + wheat_names, wheat_names=wheat_names, legume_names=[legume_name])

    planter = Planter(indexer=indexer, 
                      generation_type='forced',
                      save_wheat_positions=True,
                      positions_grid=positions_grid,
                      domain=domain)

    # Initialisation de l'instance unique de luzerne
    legume = L_egume_wrapper(
        name=legume_name,
        indexer=indexer,
        in_folder=in_folder_legume,
        out_folder=out_folder, 
        nameconfigfile='liste_usms_couplage.xls',
        ongletconfigfile='LUBBAC',
        IDusm=id_usm, 
        planter=planter, 
        caribu_scene=True
    )
       
    # Initialisation des instances multiples de blé
    wheats = []
    for w_name in wheat_names:
        w_out_folder = os.path.join(out_folder, w_name)
        os.makedirs(w_out_folder, exist_ok=True)
        
        w = Wheat_wrapper(
            name=w_name,
            in_folder=in_folder_wheat,
            out_folder=w_out_folder,
            planter=planter,
            indexer=indexer,
            external_soil_model=True,
            nitrates_uptake_forced=False,
            tillers_replications=tillers_replications,
            update_parameters_all_models=RERmax_vegetative_stages_example,
            METEO_FILENAME='LUBBAC_H_24_25_capteurs.csv',
            SENESCWHEAT_TIMESTEP=senescwheat_timestep,
            LIGHT_TIMESTEP=light_timestep,
            SOIL_PARAMETERS_FILENAME="inputs_soil_legume/Parametres_plante_exemple.xls",
            AXES_INITIAL_STATE_FILENAME="axes_initial_state.csv",
            HIDDENZONES_INITIAL_STATE_FILENAME="hiddenzones_initial_state.csv",
            rootdistribtype=wheatroot_type,
            gaic=0.2
        )
        wheats.append(w)

    lighting = Light_wrapper(
        lightmodel="caribu", 
        out_folder=out_folder, 
        sky=sky,
        planter=planter, 
        indexer=indexer,
        legume_wrapper=legume,
        writegeo=writegeo,
        geostep=geostep,
    )

    soil = Soil_wrapper(in_folder=in_folder_legume, 
                        out_folder=out_folder, 
                        IDusm=id_usm,
                        nameconfigfile='liste_usms_couplage.xls',
                        ongletconfigfile='LUBBAC', 
                        legume_wrapper=legume, 
                        planter=planter, 
                        opt_residu=0, 
                        save_results=True)

    histo_restrans_caribu = []

    ##################
    ### SIMULATION ###
    ##################
    try : 
        current_time_of_the_system = time.time()
        
        ref_w = wheats[0] # Utilisation d'un blé de référence pour les temps
        start_time = ref_w.meteo[ref_w.meteo['Date']==start_wheat].index[0]
        for w in wheats:
            w.start_time = start_time
            
        nb_iter = int(ref_w.meteo.loc[start_time, ["DOY"]].iloc[0] - legume.lsystem.DOYdeb)

        t_legume = 0

        # --------- SIMULATION LUZERNE SEULE ---------
        for t in range(nb_iter):
            legume.derive(t)

            lighting.writegeo=False
            if t%geostep == 0 :
                lighting.writegeo=True 

            scene_legume = legume.light_inputs(elements="triangles")
            lighting.run(scenes=[scene_legume], day=legume.doy(), parunit="RG")
            legume.light_results(legume.energy(), lighting)

            (
                N_content_roots_per_plant,
                roots_length_per_plant_per_soil_layer,
                plants_soil_parameters,
                plants_light_interception,
            ) = legume.soil_inputs()
            soil.run(
                legume.doy(),
                [N_content_roots_per_plant],
                [roots_length_per_plant_per_soil_layer],
                [plants_soil_parameters],
                [plants_light_interception],
            )
            legume.soil_results(soil.results, planter)

            histo_restrans_caribu.append(legume.res_trans.copy())
            legume.run()

            t_legume += 1
        
        lighting.i_vtk = lighting.i_vtk

        wheat_t_count = 0
        wheat_day_count = 0
        
        for w in wheats:
            w.GAIp =  lighting.light.totaltrianglearea / lighting.domain_area

        # <-- INITIALISATION DU POOL DE THREADS POUR LE BLÉ -->
        num_cores = os.cpu_count() or 4
        executor = concurrent.futures.ThreadPoolExecutor(max_workers=num_cores)

        # --------- CO-SIMULATION ---------
        for t_wheat in range(start_time, start_time + simulation_length, ref_w.SENESCWHEAT_TIMESTEP):
            wheat_t_count += 1
            print(f"cn_wheat timestep : {wheat_t_count}/{simulation_length}")
                
            activate_legume = ref_w.doy(t_wheat) != ref_w.next_day_next_hour(t_wheat)
            daylight = (t_wheat % light_timestep == 0) and (ref_w.PARi_next_hours(t_wheat) > 0)
            
            if daylight or activate_legume:
                if activate_legume:
                    legume.derive(t_legume)

                lighting.writegeo=False
                if t_legume%geostep == 0 :
                    lighting.writegeo=True 

                # Collecte de la lumière de tous les acteurs
                light_inputs_dict = {legume_name: legume.light_inputs(elements="triangles")}
                all_stems = []
                
                for w in wheats:
                    w_input, w_stems = w.light_inputs(planter)
                    light_inputs_dict[w.name] = w_input
                    all_stems.extend(w_stems)

                scenes = indexer.light_scenes_mgmt(light_inputs_dict)

                for w in wheats:
                    w.GAIp = lighting.light.totaltrianglearea / lighting.domain_area

                lighting.run(
                    scenes=scenes,
                    day=ref_w.doy(t_wheat),
                    hour=ref_w.hour(t_wheat),
                    parunit="RG",
                    stems=all_stems
                )
                
                if daylight:
                    for w in wheats:
                        w.light_results(energy=w.energy(t_wheat), lighting=lighting)

                if activate_legume:
                    legume.light_results(legume.energy(), lighting)
                    wheat_day_count += 1 

                    soil_inputs_dict = {legume_name: legume.soil_inputs()}

                    for w in wheats:
                        if w.rootdistribtype in ["bound", "profile"]:
                            min_depth_val = min_depth if min_depth is not None else 0.2
                            explo_rate = 0.01 # 1cm par jour en m
                            
                            w.rooting_depth = min(wheat_day_count * explo_rate + min_depth_val, soil.soil.dxyz[2][0] * len(soil.soil.dxyz[2]))
                            w.roots_bound = min(math.ceil(w.rooting_depth/soil.soil.dxyz[2][0]), len(soil.soil.dxyz[2]))

                        soil_inputs_dict[w.name] = w.soil_inputs(soil, planter, lighting)
                    
                    # Merge inputs pour le sol via l'indexer
                    (
                        N_content_roots_per_plant,
                        roots_length_per_plant_per_soil_layer,
                        plants_soil_parameters,
                        plants_light_interception,
                    ) = indexer.soil_inputs(soil_inputs_dict)
                    
                    soil.run(
                        legume.doy(),
                        N_content_roots_per_plant,
                        roots_length_per_plant_per_soil_layer,
                        plants_soil_parameters,
                        plants_light_interception
                    )
                    
                    for w in wheats:
                        w.soil_results(soil.results[4], planter=planter)
                        
                    legume.soil_results(soil.results)

                    histo_restrans_caribu.append(legume.res_trans.copy())
                    legume.run()
                    t_legume += 1

            # <-- EXÉCUTION PARALLÉLISÉE DE wheat.run(t) -->
            futures = [executor.submit(w.run, t_wheat) for w in wheats]
            concurrent.futures.wait(futures)

    finally :
        # Nettoyage et fermeture du pool de threads
        if 'executor' in locals():
            executor.shutdown(wait=True)
            
        execution_time = int(time.time() - current_time_of_the_system)
        print("\n" "Simulation run in {}".format(str(datetime.timedelta(seconds=execution_time))))
        
        legume.end()
        for w in wheats:
            w.end(run_postprocessing=run_postprocessing, run_graphs=run_graphs)
        soil.end()

        np.save(os.path.join(os.path.normpath(out_folder),"historique_restrans_caribu.npy"), np.array(histo_restrans_caribu))

if __name__ == "__main__":
    in_folder_legume = "inputs_soil_legume" 
    in_folder_wheat = "inputs_fspmwheat"
    out_folder = "outputs/WheatLuz_LUBBAC_emergence_rainbow"
    start_wheat='18/11/2024' 
    simulation_length = 1000
    id_usm=2
    writegeo = True
    geostep = 1
    run_postprocessing = False
    run_graphs = False
    
    simulation(in_folder_legume, in_folder_wheat, out_folder, 
               start_wheat, simulation_length, id_usm, min_depth=0.01,
               writegeo=writegeo, geostep=geostep, 
               run_graphs=run_graphs, run_postprocessing=run_postprocessing)