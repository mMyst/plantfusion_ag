from plantfusion.new_wheat_wrapper import Wheat_wrapper
from plantfusion.light_wrapper import Light_wrapper
from plantfusion.soil_wrapper import Soil_wrapper
from plantfusion.indexer import Indexer
from plantfusion.planter import Planter

import time
import datetime
import os
import math


def simulation(
    in_folder_legume, in_folder_wheat, out_folder,
    start_wheat, simulation_length, id_usm, 
    run_postprocessing=False, run_graphs=False, writegeo=False, geostep=1,
    wheatroot_type='profile',min_depth=None
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

    wheat_name = "wheat"
    
    indexer = Indexer(global_order=[wheat_name], wheat_names=[wheat_name])

    tillers_replications = {}
 
    sky = "turtle46"
    RERmax_vegetative_stages_example = {
        "elongwheat": {
            "RERmax": {5: 3.35e-06, 6: 2.1e-06, 7: 2.0e-06, 8: 1.83e-06, 9: 1.8e-06, 10: 1.65e-06, 11: 1.56e-06}
        }
    }
    senescwheat_timestep = 1
    light_timestep = 4
  
    ###NEW METHOD : FORCED PLANTER
    
        # Définir les paramètres d'entrée
    col_pattern = (wheat_name, "inter_row", wheat_name, wheat_name, "inter_row", wheat_name)
  
    n_rows = 6
    n_cols = 6
    cell_size = 0.05 # Taille de la cellule en mètres
    
    offset = {wheat_name: 0.15*2}
    noise = {wheat_name:0.001}

    
    planter = Planter(indexer=indexer, 
                      generation_type='grid_forced',
                      n_rows=n_rows,
                      n_cols= n_cols,
                      cell_size=cell_size,
                      col_pattern=col_pattern,
                      save_wheat_positions=True)


    wheat = Wheat_wrapper(
        in_folder=in_folder_wheat,
        out_folder=out_folder,
        planter=planter,
        indexer=indexer,
        external_soil_model=True,
        nitrates_uptake_forced=False,
        tillers_replications=tillers_replications,
        update_parameters_all_models=RERmax_vegetative_stages_example,
        METEO_FILENAME='LUBBAC_H_24_25.csv',
        SENESCWHEAT_TIMESTEP=senescwheat_timestep,
        LIGHT_TIMESTEP=light_timestep,
        SOIL_PARAMETERS_FILENAME="inputs_soil_legume/Parametres_plante_exemple.xls",
        rootdistribtype=wheatroot_type
    )

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

        if start_wheat is not None: 
            wheat.start_time=wheat.meteo[wheat.meteo['Date']==start_wheat].index[0]
        
        day_count = 0

        for t in range(wheat.start_time, wheat.start_time + simulation_length, wheat.SENESCWHEAT_TIMESTEP):



            if (bool(wheat.g.property('geometry')) and (((t % light_timestep == 0) and (wheat.PARi_next_hours(t) > 0)) or (wheat.doy(t) != wheat.next_day_next_hour(t)))):
                wheat_input, stems = wheat.light_inputs(planter)

                if  ((writegeo==True) and (t%geostep*light_timestep == 0 )) :
                    lighting.writegeo=True 
                else:
                    lighting.writegeo=False
                

                lighting.run(scenes=[wheat_input], day=wheat.doy(t), hour=wheat.hour(t), parunit="micromol.m-2.s-1", stems=stems)
                
                if ((t % light_timestep == 0) and (wheat.PARi_next_hours(t) > 0)) :
                    wheat.light_results(energy=wheat.energy(t), lighting=lighting)

                if (wheat.doy(t)  != wheat.next_day_next_hour(t) ) :
                    
                    day_count += 1

                    if wheat.rootdistribtype == "bound" or wheat.rootdistribtype == "profile":
                        min_depth = min_depth if min_depth is not None else 0.2 # unit : m
                        explo_rate = 0.01 #1cm par jour en m, approx. from Kirkegaard et Lillet 2007

                        wheat.rooting_depth = min(day_count * explo_rate + min_depth, soil.soil.dxyz[2][0] *len(soil.soil.dxyz[2]))
                        wheat.roots_bound = min(math.ceil(wheat.rooting_depth/soil.soil.dxyz[2][0]), len(soil.soil.dxyz[2])) #renvoie la couche jusqu'à laquelle les racines peuvent aller
                

                    (
                        N_content_roots_per_plant,
                        roots_length_per_plant_per_soil_layer,
                        wheat_soil_parameters,
                        plants_light_interception,
                    ) = wheat.soil_inputs(soil, planter, lighting)


                    

                    soil.run(
                        wheat.doy(t, soil3ds=True),
                        [N_content_roots_per_plant],
                        [roots_length_per_plant_per_soil_layer],
                        [wheat_soil_parameters],
                        [plants_light_interception],
                    )
                    wheat.soil_results(soil.results[4])

            wheat.run(t)

        execution_time = int(time.time() - current_time_of_the_system)
        print("\n" "Simulation run in {}".format(str(datetime.timedelta(seconds=execution_time))))

    finally:
        wheat.end(run_postprocessing=run_postprocessing, run_graphs=run_graphs)

        soil.end()


if __name__ == "__main__":
    in_folder_legume = "inputs_soil_legume"
    in_folder_wheat = "inputs_fspmwheat"
    out_folder = "outputs/wheat_LUBBAC_emergence"
    start_wheat='18/11/2024' #semis au 18/11/2024, 3 feuilles au 07/01/2025 d'après données. 31/12/2024 pour éviter pb doy, t init = 2904
    simulation_length = 4000
    id_usm = 2 #1 with reg, 2 without reg, 3 without reg and default aflalfa instead of timbale, all with perfect irrigation => only relevant for soil parameters here
    writegeo = True
    geostep = 10
    run_postprocessing = True
    run_graphs = True
    simulation(in_folder_legume, in_folder_wheat, out_folder, 
               start_wheat, simulation_length, id_usm, min_depth=0.2,
               writegeo=writegeo, geostep=geostep, 
               run_postprocessing=run_postprocessing, run_graphs=run_graphs)
