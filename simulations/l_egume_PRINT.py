from plantfusion.indexer import Indexer
from plantfusion.l_egume_wrapper import L_egume_wrapper
from plantfusion.light_wrapper import Light_wrapper
from plantfusion.soil_wrapper import Soil_wrapper
from plantfusion.planter import Planter
from plantfusion.utils import coord_from_image

from PIL import Image


from pathlib import Path
import time
import datetime
import os


def simulation(in_folder, out_folder, id_usm, writegeo, image_path):
    try:
        # Create target Directory
        os.mkdir(os.path.normpath(out_folder))
        print("Directory ", os.path.normpath(out_folder), " Created ")
    except FileExistsError:
        print("Directory ", os.path.normpath(out_folder), " already exists")

    
    legume_name = "legume"
    indexer = Indexer(global_order=[legume_name], legume_names=[legume_name])
    
    ###NEW METHOD : FORCED PLANTER
    
    image=Image.open(image_path)

    mod_img, img_coords = coord_from_image(image, n_colors=2, n_side=300,convert_mode='nb',scale=10**-3)


    # Rename the first key to 'legume'
    old_key_l = list(img_coords.keys())[1]
    old_key_w = list(img_coords.keys())[0]
    img_coords[legume_name] = img_coords.pop(old_key_l)
    #img_coords['legumeB'] = img_coords.pop(old_key_w)
    #delete any item that doesn't have a name in indexer.global_order to get rid of unused coordinates
    img_coords = {k: v for k, v in img_coords.items() if k in indexer.global_order}

    mod_img.show()
    
    planter = Planter(indexer=indexer, 
                      generation_type='forced',
                      positions_grid = img_coords)

    
    legume = L_egume_wrapper(
        name=legume_name, 
        indexer=indexer, 
        ongletconfigfile=onglet,
        in_folder=in_folder, 
        out_folder=out_folder, 
        IDusm=id_usm, 
        planter=planter,
        caribu_scene= True
    )
    
    lighting = Light_wrapper(
        lightmodel="caribu", 
        out_folder=out_folder, 
        sky='turtle46',
        planter=planter, 
        indexer=indexer,
        legume_wrapper=legume,
        writegeo=writegeo,
        infinite=False
    )

    soil = Soil_wrapper(out_folder=out_folder, ongletconfigfile=onglet, legume_wrapper=legume, legume_pattern=True)

    try:
        current_time_of_the_system = time.time()
        for t in range(legume.lsystem.derivationLength):

            legume.derive(t)

            scene_legume = legume.light_inputs(elements="triangles")
            lighting.run(scenes=[scene_legume], energy=legume.energy(), day=legume.doy(), parunit="RG")
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
    
    finally:
        execution_time = int(time.time() - current_time_of_the_system)
        print("\n" "Simulation run in {}".format(str(datetime.timedelta(seconds=execution_time))))

        legume.end()


if __name__ == "__main__":
    in_folder = "inputs_soil_legume"
    image_path = 'C:\\Users\\agrumel\\Pictures\\Saved Pictures\\newyear.png'
    out_folder = "outputs/legume_PRINT/"+ Path(image_path).stem
    onglet='LUCOS'
    id_usm = 1
    writegeo = True
    


    simulation(in_folder, out_folder, id_usm, writegeo, image_path)
