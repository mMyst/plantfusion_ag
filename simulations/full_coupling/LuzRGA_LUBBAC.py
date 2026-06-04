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



def simulation(in_folder, onglet, config_file, out_folder, id_usm, rga_usm, write_geo=False, geostep=1):
    try:
        # Create target Directory
        os.mkdir(os.path.normpath(out_folder))
        print("Directory ", os.path.normpath(out_folder), " Created ")
    except FileExistsError:
        print("Directory ", os.path.normpath(out_folder), " already exists")

    
    ###INDEXER
    legume_name = "legume"
    RGA_name = "RGA"


    indexer = Indexer(global_order=[legume_name,RGA_name], legume_names=[legume_name,RGA_name])


    # lumiere avec caribu
    sky = "turtle46"


    ###PLANTER
           # Définir les paramètres d'entrée
    # col_pattern = ("inter_row", RGA_name,legume_name, legume_name, RGA_name, "inter_row")
    col_pattern = (legume_name, RGA_name, "inter_row", "inter_row", RGA_name, legume_name)

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


    legume = L_egume_wrapper(
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
    
    rga = L_egume_wrapper(
        name=RGA_name, 
        indexer=indexer, 
        in_folder=in_folder, 
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
        writegeo=write_geo,
        geostep=geostep,
    )

    soil = Soil_wrapper(in_folder = in_folder,
                               out_folder=out_folder, 
                               nameconfigfile= config_file,
                               IDusm= id_usm,
                               ongletconfigfile= onglet,
                               planter=planter,
                               save_results= True)


    try:
        current_time_of_the_system = time.time()
        legume.lsystem.DOYdeb
        #Alfalfa headstart
        for t in range(rga.lsystem.DOYdeb - legume.lsystem.DOYdeb):
        
            legume.derive(t)

            lighting.writegeo=False
            if t%geostep == 0 :
                lighting.writegeo=True 

            
            ### CARIBU
            scene_legume = legume.light_inputs(elements="triangles")
            start = time.time()
            lighting.run(
                scenes=[scene_legume], day=legume.doy(), parunit="RG"
            )
            caribu_time = time.time() - start
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


            print("Lighting running time | CARIBU: ", caribu_time)

        #RGA 
        for t_rga in range(rga.lsystem.derivationLength):
        
            t_legume = t_rga + (rga.lsystem.DOYdeb - legume.lsystem.DOYdeb)

            legume.derive(t_legume)
            rga.derive(t_rga)

            lighting.writegeo=False
            if t%geostep == 0 :
                lighting.writegeo=True 

            
            ### CARIBU
            scene_legume = legume.light_inputs(elements="triangles")
            scene_rga = rga.light_inputs(elements="triangles")

            start = time.time()
            lighting.run(
                scenes=[scene_legume, scene_rga], day=legume.doy(), parunit="RG"
            )
            caribu_time = time.time() - start
            legume.light_results(legume.energy(), lighting)
            rga.light_results(rga.energy(), lighting)


            soil_legume_inputs = legume.soil_inputs()
            soil_rga_inputs = rga.soil_inputs()
            (
                N_content_roots_per_plant,
                roots_length_per_plant_per_soil_layer,
                plants_soil_parameters,
                plants_light_interception,
            ) = indexer.soil_inputs({legume_name : soil_legume_inputs, RGA_name : soil_rga_inputs})

            soil.run(
                legume.doy(),
                N_content_roots_per_plant,
                roots_length_per_plant_per_soil_layer,
                plants_soil_parameters,
                plants_light_interception,
            )
            legume.soil_results(soil.results, planter)
            rga.soil_results(soil.results, planter)

            legume.run()
            rga.run()


            print("Lighting running time | CARIBU: ", caribu_time)

        execution_time = int(time.time() - current_time_of_the_system)
        print("\n" "Simulation run in {}".format(str(datetime.timedelta(seconds=execution_time))))

    finally:
        legume.end()
        rga.end()
        soil.end()  


if __name__ == "__main__":
    in_folder = "inputs_soil_legume"
    out_folder = "outputs/LuzRGA_LUBBAC"
    config_file = 'liste_usms_couplage.xls'
    onglet='LUBBAC' #repiquage le 30/09, départ de la sim
    id_usm=2 #1 Timbale reg, 2 Timbale unreg
    rga_usm=7 #RGA unreg 7, RGA reg 8 
    write_geo=True
    geostep=10

    simulation(in_folder, onglet, config_file, out_folder, id_usm, rga_usm, write_geo=write_geo, geostep=geostep)
    #stade 1F le 14/10 => caler le semis en fonction


