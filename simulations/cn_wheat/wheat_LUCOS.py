from plantfusion.wheat_wrapper import Wheat_wrapper, import_fertilization
from plantfusion.light_wrapper import Light_wrapper
from plantfusion.indexer import Indexer
from plantfusion.planter import Planter


import time
import datetime
import os
import pandas as pd 


def simulation(in_folder, out_folder, start_wheat, simulation_length, mgmt, write_geo=False, run_postprocessing=False):
    try:
        # Create target Directory
        os.mkdir(os.path.normpath(out_folder))
        print("Directory ", os.path.normpath(out_folder), " Created ")
    except FileExistsError:
        print("Directory ", os.path.normpath(out_folder), " already exists")

    plants_name = "wheat"
    index_log = Indexer(global_order=[plants_name], wheat_names=[plants_name])

    N_fertilizations = import_fertilization(mgmt[0],mgmt[1],'FertNH4','inputs_fspmwheat\Lusignan_H_21-23.csv')

    #N_fertilizations = {2016: 357143, 2520: 1000000}
    


    tillers_replications = {"T1": 0.5, "T2": 0.5, "T3": 0.5, "T4": 0.5}
    plant_density = {1: 250}
    sky = [4, 5, "soc"] 
    RERmax_vegetative_stages_example = {
        "elongwheat": {
            "RERmax": {5: 3.35e-06, 6: 2.1e-06, 7: 2.0e-06, 8: 1.83e-06, 9: 1.8e-06, 10: 1.65e-06, 11: 1.56e-06}
        }
    }
    senescwheat_timestep = 1
    light_timestep = 4

    planter = Planter(generation_type="default", indexer=index_log, inter_rows=0.14, plant_density=plant_density)
    
    wheat = Wheat_wrapper(
        in_folder=in_folder,
        out_folder=out_folder,
        planter=planter,
        indexer=index_log,
        external_soil_model=False,
        nitrates_uptake_forced=False,
        N_fertilizations=N_fertilizations,
        METEO_FILENAME='Lusignan_H_21-23.csv',
        tillers_replications=tillers_replications,
        update_parameters_all_models=RERmax_vegetative_stages_example,
        SENESCWHEAT_TIMESTEP=senescwheat_timestep,
        LIGHT_TIMESTEP=light_timestep,
    )

    lighting = Light_wrapper(
        lightmodel="caribu", 
        out_folder=out_folder, 
        sky=sky,
        planter=planter, 
        indexer=index_log,
        writegeo=write_geo
    )

    current_time_of_the_system = time.time()
    wheat.start_time= wheat.meteo[wheat.meteo['Date']==start_wheat].index[0]


    for t in range(wheat.start_time, wheat.start_time+simulation_length, wheat.SENESCWHEAT_TIMESTEP):
        if (t % light_timestep == 0) and (wheat.PARi_next_hours(t) > 0):
            wheat_input, stems = wheat.light_inputs(planter)
            lighting.run(scenes=[wheat_input], day=wheat.doy(t), hour=wheat.hour(t), parunit="micromol.m-2.s-1", stems=stems)
            wheat.light_results(energy=wheat.energy(t), lighting=lighting)

        wheat.run(t)

    execution_time = int(time.time() - current_time_of_the_system)
    print("\n" "Simulation run in {}".format(str(datetime.timedelta(seconds=execution_time))))


    wheat.end(run_postprocessing=run_postprocessing)


if __name__ == "__main__":
    in_folder = "inputs_fspmwheat"
    out_folder = "outputs/cnwheat_LUCOS"
    mgmt=('inputs_soil_legume\management_LUCOS.xls','LUCOS_nomgmt') #LUCOS_highN or LUCOS_nomgmt for N+ or N- scenarios
    start_wheat= '17/12/2021'  
    simulation_length = 2500
    write_geo = True

    simulation(in_folder, out_folder, start_wheat, simulation_length, mgmt, write_geo=write_geo)
