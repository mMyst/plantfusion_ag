from plantfusion.new_wheat_wrapper import Wheat_wrapper
from plantfusion.light_wrapper import Light_wrapper
from plantfusion.soil_wrapper import Soil_wrapper
from plantfusion.indexer import Indexer
from plantfusion.planter import Planter, grid_pattern_per_cell

import time
import datetime
import os
import math
import concurrent.futures # <-- IMPORT POUR LA PARALLÉLISATION

def simulation(
    in_folder_legume, in_folder_wheat, out_folder,
    start_wheat, simulation_length, id_usm, 
    run_postprocessing=False, run_graphs=False, writegeo=False, geostep=1,
    wheatroot_type='profile',min_depth=None,
    coef_delay_til = 2.0, coef_buffer_til = 0.5, gaic = 0.11
):
    try:
        os.mkdir(os.path.normpath(out_folder))
        print("Directory ", os.path.normpath(out_folder), " Created ")
    except FileExistsError:
        print("Directory ", os.path.normpath(out_folder), " already exists")

    ######################
    ### INITIALIZATION ###
    ######################

    base_wheat_name = "wheat"
    n_rows = 6
    n_cols = 6
    cell_size = 0.05 
    
    col_pattern_base = (base_wheat_name, "inter_row", base_wheat_name, base_wheat_name, "inter_row", base_wheat_name)
    noise = {base_wheat_name: 0.001}
  
    # 1. Utilisation de la fonction externe issue de planter.py
    domain, full_grid = grid_pattern_per_cell(
        n_cols=n_cols, 
        n_rows=n_rows, 
        cell_size=cell_size, 
        col_pattern=col_pattern_base, 
        noise=noise
    )

    # 2. Filtrage pour l'Indexer : on ne garde que les clés correspondant au blé
    wheat_names = [name for name in full_grid.keys() if name.startswith(f"{base_wheat_name}_")]
    positions_grid = {name: full_grid[name] for name in wheat_names}
    
    indexer = Indexer(global_order=wheat_names, wheat_names=wheat_names)

    tillers_replications = {"T1": 0.5, "T2": 0.5, "T3": 0.5, "T4": 0.5}
    sky = "turtle46"
    RERmax_vegetative_stages_example = {
        "elongwheat": {
            "RERmax": {5: 3.35e-06, 6: 2.1e-06, 7: 2.0e-06, 8: 1.83e-06, 9: 1.8e-06, 10: 1.65e-06, 11: 1.56e-06}
        }
    }
    senescwheat_timestep = 1
    light_timestep = 4
  
    # 3. Initialisation du Planter en mode 'forced'
    planter = Planter(indexer=indexer, 
                      generation_type='forced',
                      positions_grid=positions_grid,
                      xy_plane=domain,
                      save_wheat_positions=True)

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
            rootdistribtype=wheatroot_type,
            coef_delay_til=coef_delay_til,
            coef_buffer_til=coef_buffer_til,
            gaic=gaic
        )
        wheats.append(w)

    lighting = Light_wrapper(
        lightmodel="caribu", 
        out_folder=out_folder, 
        sky=sky,
        planter=planter, 
        indexer=indexer,
        writegeo=writegeo
    )

    soil = Soil_wrapper(in_folder=in_folder_legume, 
                        out_folder=out_folder, 
                        nameconfigfile='liste_usms_couplage.xls',
                        IDusm=id_usm,
                        ongletconfigfile='LUBBAC',
                        planter=planter, 
                        opt_residu=0, 
                        save_results=True)
    
    ##################
    ### SIMULATION ###
    ##################

    try :
        current_time_of_the_system = time.time()
        
        ref_w = wheats[0] 

        if start_wheat is not None:
            start_time = ref_w.meteo[ref_w.meteo['Date']==start_wheat].index[0]
            for w in wheats:
                w.start_time = start_time
        else:
            start_time = ref_w.start_time
        
        day_count = 0
        
        # <-- INITIALISATION DU POOL DE THREADS -->
        num_cores = os.cpu_count() or 4
        executor = concurrent.futures.ThreadPoolExecutor(max_workers=num_cores)

        for t in range(start_time, start_time + simulation_length, ref_w.SENESCWHEAT_TIMESTEP):

            if (bool(ref_w.g.property('geometry')) and (((t % light_timestep == 0) and (ref_w.PARi_next_hours(t) > 0)) or (ref_w.doy(t) != ref_w.next_day_next_hour(t)))):
                
                light_inputs_dict = {}
                all_stems = []
                for w in wheats:
                    w_input, w_stems = w.light_inputs(planter)
                    light_inputs_dict[w.name] = w_input
                    all_stems.extend(w_stems)

                if ((writegeo==True) and (t%geostep*light_timestep == 0 )) :
                    lighting.writegeo=True 
                else:
                    lighting.writegeo=False
                
                merged_scene = indexer.light_scenes_mgmt(light_inputs_dict)
                lighting.run(scenes=merged_scene, day=ref_w.doy(t), hour=ref_w.hour(t), parunit="micromol.m-2.s-1", stems=all_stems)
                
                if ((t % light_timestep == 0) and (ref_w.PARi_next_hours(t) > 0)) :
                    for w in wheats:
                        w.light_results(energy=w.energy(t), lighting=lighting)

                if (ref_w.doy(t) != ref_w.next_day_next_hour(t)) :
                    
                    day_count += 1
                    
                    N_content_roots_all = []
                    roots_length_all = []
                    wheat_soil_params_all = []
                    plants_light_interception_all = []

                    for w in wheats:
                        if w.rootdistribtype == "bound" or w.rootdistribtype == "profile":
                            min_depth_val = min_depth if min_depth is not None else 0.2
                            explo_rate = 0.01

                            w.rooting_depth = min(day_count * explo_rate + min_depth_val, soil.soil.dxyz[2][0] *len(soil.soil.dxyz[2]))
                            w.roots_bound = min(math.ceil(w.rooting_depth/soil.soil.dxyz[2][0]), len(soil.soil.dxyz[2]))
                    
                        (
                            N_content_roots,
                            roots_length,
                            wheat_soil_parameters,
                            plants_light_interception,
                        ) = w.soil_inputs(soil, planter, lighting)
                        
                        N_content_roots_all.append(N_content_roots)
                        roots_length_all.append(roots_length)
                        wheat_soil_params_all.append(wheat_soil_parameters)
                        plants_light_interception_all.append(plants_light_interception)

                    soil.run(
                        ref_w.doy(t, soil3ds=True),
                        N_content_roots_all,
                        roots_length_all,
                        wheat_soil_params_all,
                        plants_light_interception_all,
                    )
                    
                    for w in wheats:
                        w.soil_results(soil.results[4])

            print(f"Step {t-start_time}/{simulation_length}")
            
            # <-- EXÉCUTION PARALLÉLISÉE DE wheat.run(t) -->
            futures = [executor.submit(w.run, t) for w in wheats]
            concurrent.futures.wait(futures)
            
    finally:
        # Nettoyage et fermeture du pool de threads
        if 'executor' in locals():
            executor.shutdown(wait=True)
            
        execution_time = int(time.time() - current_time_of_the_system)
        print("\n" "Simulation run in {}".format(str(datetime.timedelta(seconds=execution_time))))
        
        print(wheats[0].tillers_replications)
        soil.end()
        for w in wheats:
            w.end(run_postprocessing=run_postprocessing, run_graphs=run_graphs)
        
if __name__ == "__main__":
    in_folder_legume = "inputs_soil_legume"
    in_folder_wheat = "inputs_fspmwheat"
    out_folder = "outputs/wheat_rainbow_parallel"
    start_wheat='18/11/2024' 
    simulation_length = 500 
    id_usm = 2 
    writegeo = True
    geostep = 1
    run_postprocessing = True
    run_graphs = True
    simulation(in_folder_legume, in_folder_wheat, out_folder, 
               start_wheat, simulation_length, id_usm, min_depth=0.01,
               writegeo=writegeo, geostep=geostep, 
               coef_delay_til = 2.0, coef_buffer_til = 0.5, gaic = 0.20,
               run_postprocessing=run_postprocessing, run_graphs=run_graphs)