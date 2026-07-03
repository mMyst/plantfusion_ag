from plantfusion.l_egume_wrapper import L_egume_wrapper
from plantfusion.light_wrapper import Light_wrapper
from plantfusion.soil_wrapper import Soil_wrapper
from plantfusion.planter import Planter
from plantfusion.indexer import Indexer
from plantfusion.utils import create_child_folder

import numpy as np


import os
import time
import datetime
import pandas


def simulation(in_folder, out_folder, id_usm, write_geo=False):
    try:
        # Create target Directory
        os.mkdir(os.path.normpath(out_folder))
        print("Directory ", os.path.normpath(out_folder), " Created ")
    except FileExistsError:
        print("Directory ", os.path.normpath(out_folder), " already exists")


    plants_name = "legume"
    index_log = Indexer(global_order=[plants_name], legume_names=[plants_name])

    #sky = "turtle46"
    sky = "inputs_soil_legume/sky_5.data"

    # version par défaut
    planter = Planter(indexer=index_log, legume_cote={plants_name : 40.}, legume_number_of_plants={plants_name : 64})

    # lumiere avec caribu
    legume_caribu = L_egume_wrapper(
        name=plants_name, 
        indexer=index_log, 
        in_folder=in_folder, 
        out_folder= out_folder, 
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
    )
    soil_caribu = Soil_wrapper(out_folder=out_folder, legume_wrapper=legume_caribu,  legume_pattern=True, planter=planter)

    light_data = {"epsi": [], "parip": [], "t": []}

    histo_restrans_caribu = []

    try:
        current_time_of_the_system = time.time()
        for t in range(legume_caribu.lsystem.derivationLength):
            legume_caribu.derive(t)


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

            newpars_caribu = legume_caribu.res_trans[-1]/legume_caribu.lsystem.tag_loop_inputs[15]
            
            histo_restrans_caribu.append(legume_caribu.res_trans.copy())
            legume_caribu.run()

            print("Lighting running time |  CARIBU: ", caribu_time)

        execution_time = int(time.time() - current_time_of_the_system)
        print("\n" "Simulation run in {}".format(str(datetime.timedelta(seconds=execution_time))))

    finally:
        legume_caribu.end()
        np.save(os.path.join(os.path.normpath(out_folder),"historique_restrans_caribu.npy"), np.array(histo_restrans_caribu))



if __name__ == "__main__":
    in_folder = "inputs_soil_legume"
    out_folder = "outputs/legume_caribu_buggrid"
    write_geo = False

    simulation(in_folder, out_folder, 1711, write_geo)
