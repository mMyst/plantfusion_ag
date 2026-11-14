from plantfusion.adel_wheat_wrapper import Wheat_wrapper
from plantfusion.l_egume_wrapper import L_egume_wrapper
from plantfusion.light_wrapper import Light_wrapper
from plantfusion.soil_wrapper import Soil_wrapper
from plantfusion.indexer import Indexer
from plantfusion.planter import Planter, grid_pattern
from plantfusion.utils import get_top_height
from elongwheat import tillering as elongwheat_tillering

import numpy as np
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
    legume_name = "legume"
    indexer = Indexer(global_order=[legume_name, wheat_name], wheat_names=[wheat_name], legume_names=[legume_name])

    tillers_replications = {"T1": 0.5, "T2": 0.5, "T3": 0.5, "T4": 0.5, "T5":0.5}
 
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
    col_pattern = (legume_name, "inter_row", wheat_name, wheat_name, "inter_row", legume_name)
  
    n_rows = 6
    n_cols = 6
    cell_size = 0.05 # Taille de la cellule en mètres
    
    offset = {wheat_name: 0.15*2, legume_name: 0}
    noise = {wheat_name:0.001,legume_name:0.01}


    domain, full_grid = grid_pattern(
        n_cols=n_cols, 
        n_rows=n_rows, 
        cell_size=cell_size, 
        col_pattern=col_pattern, 
        noise=None
    )

    
    planter = Planter(indexer=indexer, 
                      generation_type='forced',
                      save_wheat_positions=True,
                      positions_grid=full_grid,
                      domain = domain)


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

       
    wheat = Wheat_wrapper(
        in_folder=in_folder_wheat,
        out_folder=out_folder,
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
        gaic=0.07
    )

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
        t_legume = 0
        wheat.start_time=wheat.meteo[wheat.meteo['Date']==start_wheat].index[0]
        nb_iter = int(wheat.meteo.loc[wheat.start_time, ["DOY"]].iloc[0] - legume.lsystem.DOYdeb)

        # onlylegume_index = 0
        # save_planter_nb_plants = planter.number_of_plants
        # planter.number_of_plants[onlylegume_index] = planter.number_of_plants[legume.global_index]
        
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

        # planter.number_of_plants = save_planter_nb_plants
        
        lighting.i_vtk = lighting.i_vtk

        wheat_t_count = 0
        wheat_day_count = 0
        

        #compute GAIp using total triangle area stored in LVM and domain area from planter
        # wheat.GAIp =  lighting.light.totaltrianglearea / lighting.domain_area
        # wheat.GAIp = elongwheat_tillering.compute_local_GAIp(
        #     scene = lighting.light._LightVegeManager__complete_trimesh, #__geometry["scenes"],
        #     radius = 0.10,
        #     target_point=(0.15, 0.15, 0.0),
        #     xy_only = True
        # )

        for t_wheat in range(wheat.start_time,wheat.start_time +simulation_length, wheat.SENESCWHEAT_TIMESTEP):
            wheat_t_count += 1
            print("cn_wheat timestep : "+str(wheat_t_count)+"/"+str(simulation_length))

                
            activate_legume = wheat.doy(t_wheat) != wheat.next_day_next_hour(t_wheat)
            daylight = (t_wheat % light_timestep == 0) and (wheat.PARi_next_hours(t_wheat) > 0)
            

            if daylight or activate_legume:
           
                
                if activate_legume:
                    legume.derive(t_legume)

                lighting.writegeo=False
                if t_legume%geostep == 0 :
                    lighting.writegeo=True 

                wheat_input, stems = wheat.light_inputs(planter)
                legume_input = legume.light_inputs(elements="triangles")
                scenes = indexer.light_scenes_mgmt({wheat_name : wheat_input, legume_name : legume_input})

                #check legume height vs wheat height and schedule a cut if needed
                #/!\ legume needs to be triangles 
                w_height = get_top_height(wheat_input)
                l_height = get_top_height(legume_input)/100




                # wheat.GAIp = lighting.light.totaltrianglearea / lighting.domain_area
                wheat.GAIp = elongwheat_tillering.compute_local_GAIp(
                    scene = lighting.light.to_plantGL(),
                    radius = 0.10,
                    target_point=(0.125, 0.125, 0.0),
                    xy_only = True
                )

                test_gai = (wheat.GAIp > wheat.GAIc-0.005 and wheat_t_count < 1000)
                test_height = (l_height > 0.10 and l_height > w_height)

                if test_gai or test_height : #l_height > 0.10 and l_height > w_height:
                    print("✂️ LET'S GET THE SHEARS OUT ✂️")
                    legume.schedule_a_cut(1)

                lighting.run(
                    scenes=scenes,
                    day=wheat.doy(t_wheat),
                    hour=wheat.hour(t_wheat),
                    parunit="RG",
                    stems=stems
                )
                if daylight:
                    wheat.light_results(energy=wheat.energy(t_wheat), lighting=lighting)

                if activate_legume:
                    legume.light_results(legume.energy(), lighting)

                    wheat_day_count += 1 

                    if wheat.rootdistribtype == "bound" or wheat.rootdistribtype == "profile":
                        min_depth = min_depth if min_depth is not None else 0.2 # unit : m
                        explo_rate = 0.01 #1cm par jour en m, approx. from Kirkegaard et Lillet 2007

                        wheat.rooting_depth = min(wheat_day_count * explo_rate + min_depth, soil.soil.dxyz[2][0] *len(soil.soil.dxyz[2]))
                        wheat.roots_bound = min(math.ceil(wheat.rooting_depth/soil.soil.dxyz[2][0]), len(soil.soil.dxyz[2])) #renvoie la couche jusqu'à laquelle les racines peuvent aller
        


                    soil_wheat_inputs = wheat.soil_inputs(soil, planter, lighting)
                    soil_legume_inputs = legume.soil_inputs()
                    (
                        N_content_roots_per_plant,
                        roots_length_per_plant_per_soil_layer,
                        plants_soil_parameters,
                        plants_light_interception,
                    ) = indexer.soil_inputs({legume_name : soil_legume_inputs, wheat_name : soil_wheat_inputs})
                    
                    soil.run(
                        legume.doy(),
                        N_content_roots_per_plant,
                        roots_length_per_plant_per_soil_layer,
                        plants_soil_parameters,
                        plants_light_interception
                    )
                    wheat.soil_results(soil.results[4], planter=planter)
                    legume.soil_results(soil.results)

                    histo_restrans_caribu.append(legume.res_trans.copy())
                    legume.run()
                    t_legume += 1

            wheat.run(t_wheat)

        execution_time = int(time.time() - current_time_of_the_system)
        print("\n" "Simulation run in {}".format(str(datetime.timedelta(seconds=execution_time))))

    finally :
        legume.end()
        wheat.end(run_postprocessing=run_postprocessing, run_graphs=run_graphs)
        soil.end()

        np.save(os.path.join(os.path.normpath(out_folder),"historique_restrans_caribu.npy"), np.array(histo_restrans_caribu))



if __name__ == "__main__":
    in_folder_legume = "inputs_soil_legume" #debut luzerne le 30/09/2024 => avancé au 24 septembre (doy 269) pour fitter l'apparition des feuilles
    in_folder_wheat = "inputs_fspmwheat"
    out_folder = "outputs/WheatLuzREG_LUBBAC_emergence_adel_autoGAI"
    start_wheat='18/11/2024' #semis au 18/11/2024, 3 feuilles au 07/01/2025 d'après données.
    simulation_length = 1500 #3370 = meme cumul thermique que référence. 3430 = arrêt à la fermeture du couvert blé d'après scan eric (09/04/2025)
    id_usm=2 #1 with reg, 2 without reg, 3 without reg and default aflalfa instead of timbale, all with perfect irrigation
    #11 reg meteo capteurs, #12 noreg  TODO : add accurate irrigation
    writegeo = True
    geostep = 10
    run_postprocessing = True
    run_graphs = True
    simulation(in_folder_legume, in_folder_wheat, out_folder, 
               start_wheat, simulation_length, id_usm, min_depth=0.01,
               writegeo=writegeo, geostep=geostep, 
               run_graphs=run_graphs, run_postprocessing=run_postprocessing)
