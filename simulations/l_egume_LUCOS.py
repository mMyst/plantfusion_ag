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



def simulation(in_folder, onglet, out_folder, id_usm, write_geo=False, geostep=1):
    try:
        # Create target Directory
        os.mkdir(os.path.normpath(out_folder))
        print("Directory ", os.path.normpath(out_folder), " Created ")
    except FileExistsError:
        print("Directory ", os.path.normpath(out_folder), " already exists")

    
    legume_name = "legume"
    sky = "turtle46"
    #"inputs_soil_legume/sky_5.data"

    index_log = Indexer(global_order=[legume_name], legume_names=[legume_name])


    # lumiere avec caribu

    ###
           # Définir les paramètres d'entrée
    densities = {legume_name: 300}
    n_rows = {legume_name: 2}
    inter_rows = {legume_name: 0.14}
    offset = {legume_name: 0}
    noise = {legume_name:0.01}

    
    planter = Planter(indexer=index_log, 
                      generation_type='row_forced',
                      plant_density=densities,
                      n_rows=n_rows,
                      inter_rows=inter_rows)


    legume_caribu = L_egume_wrapper(
        name=legume_name, 
        indexer=index_log, 
        in_folder=in_folder, 
        ongletconfigfile=onglet,
        out_folder=out_folder, 
        IDusm=id_usm, 
        caribu_scene=True, 
        planter=planter
    
    )
    
    lighting_caribu = Light_wrapper(
        lightmodel="caribu",
        out_folder=out_folder,
        indexer=index_log, 
        planter=planter, 
        legume_wrapper=legume_caribu,
        sky=sky,
        writegeo=write_geo,
        geostep=geostep,
    )

    soil_caribu = Soil_wrapper(out_folder=out_folder, legume_wrapper=legume_caribu,  legume_pattern=True, planter=planter)



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


if __name__ == "__main__":
    in_folder = "inputs_soil_legume"
    out_folder = "outputs/legume_LUCOS"
    id_usm=7
    write_geo=True
    onglet='LUCOS'
    geostep=50

    simulation(in_folder, onglet, out_folder, id_usm, write_geo=write_geo, geostep=geostep)

