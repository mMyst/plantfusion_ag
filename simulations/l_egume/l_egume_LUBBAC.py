from plantfusion.l_egume_wrapper import L_egume_wrapper, passive_lighting
from plantfusion.light_wrapper import Light_wrapper
from plantfusion.soil_wrapper import Soil_wrapper
from plantfusion.planter import Planter
from plantfusion.indexer import Indexer
from plantfusion.utils import create_child_folder

import numpy
import os
import time
import datetime



def simulation(in_folder, onglet, config_file, out_folder, id_usm, write_geo=False, geostep=1):
    try:
        # Create target Directory
        os.mkdir(os.path.normpath(out_folder))
        print("Directory ", os.path.normpath(out_folder), " Created ")
    except FileExistsError:
        print("Directory ", os.path.normpath(out_folder), " already exists")

    
    ###INDEXER
    legume_name = "legume"

    indexer = Indexer(global_order=[legume_name], legume_names=[legume_name])


    # lumiere avec caribu
    sky = "turtle46"
    sky = "inputs_soil_legume/sky_5.data"


    ###PLANTER
           # Définir les paramètres d'entrée
    col_pattern = ("inter_row", "inter_row",legume_name, legume_name, "inter_row", "inter_row")
    col_pattern = (legume_name,legume_name,legume_name,legume_name,legume_name,legume_name)

    n_rows = 6
    n_cols = 6
    cell_size = 0.05 # Taille de la cellule en mètres

    offset = {legume_name: 0}
    noise = {legume_name:0.01}

    
    planter = Planter(indexer=indexer, 
                      generation_type='grid_forced',
                      n_rows=n_rows,
                      n_cols= n_cols,
                      cell_size=cell_size,
                      col_pattern=col_pattern)


    legume_caribu = L_egume_wrapper(
        name=legume_name, 
        indexer=indexer, 
        in_folder=in_folder, 
        nameconfigfile= config_file,
        ongletconfigfile=onglet,
        out_folder=out_folder, 
        IDusm=id_usm, 
        caribu_scene=True, 
        planter=planter
    
    )
    
    lighting_caribu = Light_wrapper(
        lightmodel="caribu",
        out_folder=out_folder,
        indexer=indexer, 
        planter=planter, 
        legume_wrapper=legume_caribu,
        sky=sky,
        writegeo=write_geo,
        geostep=geostep,
    )

    soil_caribu = Soil_wrapper(in_folder = in_folder,
                               out_folder=out_folder, 
                               nameconfigfile= config_file,
                               ongletconfigfile= onglet,
                               legume_wrapper=legume_caribu,  
                               planter=planter,
                               save_results= True)



    try:
        current_time_of_the_system = time.time()
        for t in range(legume_caribu.lsystem.derivationLength):
        
            legume_caribu.derive(t)

            lighting_caribu.writegeo=False
            if t%geostep == 0 :
                lighting_caribu.writegeo=True 

            
            ### CARIBU
            scene_legume = legume_caribu.light_inputs(elements="triangles")
            start = time.time()
            lighting_caribu.run(
                scenes=[scene_legume], day=legume_caribu.doy(), parunit="RG"
            )
            caribu_time = time.time() - start
            legume_caribu.light_results(legume_caribu.energy(), lighting_caribu)

            (
                N_content_roots_per_plant,
                roots_length_per_plant_per_soil_layer,
                plants_soil_parameters,
                plants_light_interception,
            ) = legume_caribu.soil_inputs()

            soil_caribu.run(
                legume_caribu.doy(),
                [N_content_roots_per_plant],
                [roots_length_per_plant_per_soil_layer],
                [plants_soil_parameters],
                [plants_light_interception],
            )
            legume_caribu.soil_results(soil_caribu.results, planter)

            legume_caribu.run()

            print("Lighting running time | CARIBU: ", caribu_time)

        execution_time = int(time.time() - current_time_of_the_system)
        print("\n" "Simulation run in {}".format(str(datetime.timedelta(seconds=execution_time))))

    finally:
        legume_caribu.end()
        soil_caribu.end()  




if __name__ == "__main__":
    in_folder = "inputs_soil_legume"
    out_folder = "outputs/legume_LUBBAC_allfixed"
    config_file = 'liste_usms_couplage.xls'
    onglet='LUBBAC' #repiquage le 30/09, départ de la sim
    id_usm=12 #1 with reg, 2 without reg, 3 without reg and default aflalfa instead of timbale, all with perfect irrigation
    write_geo=True
    geostep=1

    simulation(in_folder, onglet, config_file, out_folder, id_usm, write_geo=write_geo, geostep=geostep)
    #stade 1F le 14/10 => caler le semis en fonction


