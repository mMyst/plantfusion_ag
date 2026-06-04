from plantfusion.wheat_wrapper import Wheat_wrapper
from plantfusion.l_egume_wrapper import L_egume_wrapper, passive_lighting
from plantfusion.light_wrapper import Light_wrapper
from plantfusion.soil_wrapper import Soil_wrapper
from plantfusion.planter import Planter
from plantfusion.indexer import Indexer

import os
import time
import datetime
import math




def simulation(
    in_folder_legume, in_folder_wheat, out_folder,
    start_wheat, simulation_length, id_usm, rga_usm,
    onglet, config_file, writegeo=False, geostep=1, 
    run_postprocessing=False, run_graphs=False,
    wheatroot_type='profile',min_depth=None
):
    try:
        # Create target Directory
        os.mkdir(os.path.normpath(out_folder))
        print("Directory ", os.path.normpath(out_folder), " Created ")
    except FileExistsError:
        print("Directory ", os.path.normpath(out_folder), " already exists")

    
    ###INDEXER
    wheat_name = "wheat"
    legume_name = "legume"
    RGA_name = "RGA"


    indexer = Indexer(global_order=[legume_name,RGA_name,wheat_name], wheat_names=[wheat_name],legume_names=[legume_name,RGA_name])

    tillers_replications = {"T1": 0.5, "T2": 0.5}

    # lumiere avec caribu
    sky = "turtle46"
    RERmax_vegetative_stages_example = {
        "elongwheat": {
            "RERmax": {5: 3.35e-06, 6: 2.1e-06, 7: 2.0e-06, 8: 1.83e-06, 9: 1.8e-06, 10: 1.65e-06, 11: 1.56e-06}
        }
    }
    senescwheat_timestep = 1
    light_timestep = 4

    ###PLANTER
           # Définir les paramètres d'entrée
    # col_pattern = (wheat_name, RGA_name,legume_name, legume_name, RGA_name, wheat_name)
    col_pattern = (legume_name, RGA_name, wheat_name, wheat_name, RGA_name, legume_name)

    n_rows = 6
    n_cols = 6
    cell_size = 0.05 # Taille de la cellule en mètres

    offset = {wheat_name: 0, legume_name: 0, RGA_name: 0}
    noise = {wheat_name:0.001,legume_name:0.01, RGA_name:0.01}


    
    planter = Planter(indexer=indexer, 
                      generation_type='grid_forced',
                      n_rows=n_rows,
                      n_cols= n_cols,
                      cell_size=cell_size,
                      col_pattern=col_pattern)


    legume = L_egume_wrapper(
        name=legume_name, 
        indexer=indexer, 
        in_folder=in_folder_legume, 
        nameconfigfile= config_file,
        ongletconfigfile=onglet,
        out_folder=out_folder, 
        IDusm=id_usm, 
        caribu_scene=True, 
        planter=planter
    
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
        METEO_FILENAME='LUBBAC_H_24_25.csv',
        SENESCWHEAT_TIMESTEP=senescwheat_timestep,
        LIGHT_TIMESTEP=light_timestep,
        SOIL_PARAMETERS_FILENAME="inputs_soil_legume/Parametres_plante_exemple.xls",
        AXES_INITIAL_STATE_FILENAME="axes_initial_state_2til.csv",
        HIDDENZONES_INITIAL_STATE_FILENAME="hiddenzones_initial_state_2til.csv",
        rootdistribtype=wheatroot_type
    )
    
    rga = L_egume_wrapper(
        name=RGA_name, 
        indexer=indexer, 
        in_folder=in_folder_legume, 
        nameconfigfile= config_file,
        ongletconfigfile=onglet,    
        out_folder=out_folder, 
        IDusm=rga_usm, 
        caribu_scene=True, 
        planter=planter
    )



    lighting = Light_wrapper(
        lightmodel="caribu",
        out_folder=out_folder,
        indexer=indexer, 
        planter=planter, 
        legume_wrapper=[legume, rga],
        sky=sky,
        writegeo=writegeo,
        geostep=geostep,
    )

    soil = Soil_wrapper(in_folder = in_folder_legume,
                               out_folder=out_folder, 
                               nameconfigfile= config_file,
                               IDusm= id_usm,
                               ongletconfigfile= onglet,
                               planter=planter,
                               save_results= True)



    ##################
    ### SIMULATION ###
    ##################

    current_time_of_the_system = time.time()
    t_legume = 0
    wheat.start_time=wheat.meteo[wheat.meteo['Date']==start_wheat].index[0]
    nb_iter = int(wheat.meteo.loc[wheat.start_time, ["DOY"]].iloc[0] - legume.lsystem.DOYdeb)

    # onlylegume_index = 0
    # save_planter_nb_plants = planter.number_of_plants
    # planter.number_of_plants[onlylegume_index] = planter.number_of_plants[legume.global_index]
    
    #only legume for alfalfa
    for t in range(nb_iter):

        legume.derive(t)

        lighting.writegeo=False
        if t%geostep == 0 :
            lighting.writegeo=True 


        lighting.writegeo=False
        if t%geostep == 0 :
            lighting.writegeo=True 

        
        ### CARIBU
        scene_legume = legume.light_inputs(elements="triangles")
        lighting.run(
            scenes=[scene_legume], day=legume.doy(), parunit="RG"
        )
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

        legume.run()

        t_legume += 1

    
    lighting.i_vtk = lighting.i_vtk

    wheat_t_count = 0
    wheat_day_count = 0

    #wheat and rga start at the same time
    for t_wheat in range(wheat.start_time,wheat.start_time +simulation_length, wheat.SENESCWHEAT_TIMESTEP):
        wheat_t_count += 1
        print("cn_wheat timestep : "+str(wheat_t_count))
              
        activate_legume = wheat.doy(t_wheat) != wheat.next_day_next_hour(t_wheat)
        daylight = (t_wheat % light_timestep == 0) and (wheat.PARi_next_hours(t_wheat) > 0)

        if daylight or activate_legume:
            if activate_legume:
                legume.derive(t_legume)
                rga.derive(t_legume)


            lighting.writegeo=False
            if t_legume%geostep == 0 :
                lighting.writegeo=True 

            wheat_input, stems = wheat.light_inputs(planter)
            legume_input = legume.light_inputs(elements="triangles")
            rga_input = rga.light_inputs(elements="triangles")
            scenes = indexer.light_scenes_mgmt({wheat_name : wheat_input, legume_name : legume_input, RGA_name : rga_input})

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
                rga.light_results(rga.energy(), lighting)

                wheat_day_count += 1 

                if wheat.rootdistribtype == "bound" or wheat.rootdistribtype == "profile":
                    min_depth = min_depth if min_depth is not None else 0.2 # unit : m
                    explo_rate = 0.01 #1cm par jour en m, approx. from Kirkegaard et Lillet 2007

                    wheat.rooting_depth = min(wheat_day_count * explo_rate + min_depth, soil.soil.dxyz[2][0] *len(soil.soil.dxyz[2]))
                    wheat.roots_bound = min(math.ceil(wheat.rooting_depth/soil.soil.dxyz[2][0]), len(soil.soil.dxyz[2])) #renvoie la couche jusqu'à laquelle les racines peuvent aller
      


                soil_wheat_inputs = wheat.soil_inputs(soil, planter, lighting)
                soil_legume_inputs = legume.soil_inputs()
                soil_rga_inputs = rga.soil_inputs()
                (
                    N_content_roots_per_plant,
                    roots_length_per_plant_per_soil_layer,
                    plants_soil_parameters,
                    plants_light_interception,
                ) = indexer.soil_inputs({legume_name : soil_legume_inputs, wheat_name : soil_wheat_inputs, RGA_name : soil_rga_inputs})
                
                soil.run(
                    legume.doy(),
                    N_content_roots_per_plant,
                    roots_length_per_plant_per_soil_layer,
                    plants_soil_parameters,
                    plants_light_interception
                )
                wheat.soil_results(soil.results[4], planter=planter)
                legume.soil_results(soil.results)
                rga.soil_results(soil.results)


                legume.run()
                rga.run()

                t_legume += 1

        wheat.run(t_wheat)

    execution_time = int(time.time() - current_time_of_the_system)
    print("\n" "Simulation run in {}".format(str(datetime.timedelta(seconds=execution_time))))


    legume.end()
    rga.end()
    wheat.end(run_postprocessing=run_postprocessing, run_graphs=run_graphs)
    soil.end()



if __name__ == "__main__":
    in_folder_legume = "inputs_soil_legume"
    in_folder_wheat = "inputs_fspmwheat/forced_tillers_init"
    out_folder = "outputs/WheatLuzRGA_LUBBAC"
    start_wheat='31/12/2024' #semis au 18/11/2024, 3 feuilles au 07/01/2025 d'après données. 31/12/2024 pour éviter pb doy, t init = 2904
    config_file = 'liste_usms_couplage.xls'
    simulation_length = 2500
    onglet='LUBBAC' #repiquage le 30/09, départ de la sim
    id_usm=2 #1 T1mbale reg, 2 Timbale unreg
    rga_usm=7 #RGA unreg 7, RGA reg 8 
    writegeo=True
    geostep=10
    run_postprocessing = True
    run_graphs = True

    simulation(in_folder_legume, in_folder_wheat, out_folder, 
               start_wheat, simulation_length, 
               id_usm, rga_usm, 
               onglet, config_file,
               wheatroot_type='profile',min_depth=0.2,
               writegeo=writegeo, geostep=geostep, 
               run_graphs=run_graphs, run_postprocessing=run_postprocessing)
    #stade 1F le 14/10 => caler le semis en fonction


