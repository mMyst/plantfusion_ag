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



def simulation(in_folder, onglet, out_folder, id_usm, write_geo=False):
    try:
        # Create target Directory
        os.mkdir(os.path.normpath(out_folder))
        print("Directory ", os.path.normpath(out_folder), " Created ")
    except FileExistsError:
        print("Directory ", os.path.normpath(out_folder), " already exists")

   
    plants_name = "legume"
    index_log = Indexer(global_order=[plants_name], legume_names=[plants_name])


    # lumiere avec caribu

    ##planter using old legume puppeting method 
    #planter = Planter(indexer=index_log, legume_cote={plants_name : 40.}, legume_number_of_plants={plants_name : 64},inter_rows=0.14)
    
    #new planter method under development 
    grid = {"legume": [[numpy.array([1.,2.,0.]), numpy.array([11.,2.,0.]), numpy.array([7.,2.,0.])]],
                "wheat":[[numpy.array([1.,7.,0.]), numpy.array([11.,7.,0.]), numpy.array([7.,7.,0.])]],
                "other":[[numpy.array([1.,12.,0.]), numpy.array([11.,12.,0.]), numpy.array([7.,12.,0.])]]
                }
    

    ###
        # Définir les paramètres d'entrée
    densities = {'blé': 150, 'luzerne': 1000,'orge':200}
    n_rows = {'blé': 4, 'luzerne': 4,'orge':4}
    inter_rows = {'blé': 0.15, 'luzerne': 0.15,'orge':0.15}
    offset = {'blé': 0.15*2, 'luzerne': 0}
    noise = {'blé':0.001,'luzerne':0.01,'orge':0.001}

    #planter = Planter(indexer=index_log, generation_type='forced',positions_grid=grid)
    planter = Planter(indexer=index_log, generation_type='row_forced',plant_density=densities)


    legume_caribu = L_egume_wrapper(
        name=plants_name, 
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
        sky="inputs_soil_legume/sky_5.data",
        writegeo=write_geo,
    )
    soil_caribu = Soil_wrapper(out_folder=out_folder, legume_wrapper=legume_caribu,  legume_pattern=True, planter=planter)



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

            legume_caribu.run()

            print("Lighting running time | CARIBU: ", caribu_time)

        execution_time = int(time.time() - current_time_of_the_system)
        print("\n" "Simulation run in {}".format(str(datetime.timedelta(seconds=execution_time))))

    finally:
        legume_caribu.end()


if __name__ == "__main__":
    in_folder = "inputs_soil_legume"
    out_folder = "outputs/legume_LUCOS"
    id_usm=6
    write_geo=True
    onglet='LUCOS'

    simulation(in_folder, onglet, out_folder, id_usm, write_geo=write_geo)
