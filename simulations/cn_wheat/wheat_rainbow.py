from plantfusion.wheat_wrapper import Wheat_wrapper
from plantfusion.light_wrapper import Light_wrapper
from plantfusion.indexer import Indexer
from plantfusion.planter import Planter

import time
import datetime
import os


def simulation(in_folder, out_folder,start_wheat=None, simulation_length=2500, write_geo=False, run_postprocessing=False, run_graphs=False):
    try:
        # Create target Directory
        os.mkdir(os.path.normpath(out_folder))
        print("Directory ", os.path.normpath(out_folder), " Created ")
    except FileExistsError:
        print("Directory ", os.path.normpath(out_folder), " already exists")

    plants_name = "wheat"
    index_log = Indexer(global_order=[plants_name], wheat_names=[plants_name])

    N_fertilizations = {2016: 357143, 2520: 1000000}
    tillers_replications = {"T1": 0.5, "T2": 0.5, "T3": 0.5, "T4": 0.5}
    plant_density = {1: 250}
    sky = "turtle46"
    RERmax_vegetative_stages_example = {
        "elongwheat": {
            "RERmax": {5: 3.35e-06, 6: 2.1e-06, 7: 2.0e-06, 8: 1.83e-06, 9: 1.8e-06, 10: 1.65e-06, 11: 1.56e-06}
        }
    }
    senescwheat_timestep = 1
    light_timestep = 4

    positions_greed = {}
    cell_size =0.5

    planter = Planter(generation_type="forced", indexer=index_log, inter_rows=0.15, plant_density=plant_density)
    wheat_dict= {}
    n_wheat = 64

    for i in range (0,n_wheat):
        wheat_dict['wheat'+(i+1)]=Wheat_wrapper(
            in_folder=in_folder,
            out_folder=out_folder,
            planter=planter,
            indexer=index_log,
            external_soil_model=False,
            nitrates_uptake_forced=False,
            N_fertilizations=N_fertilizations,
            tillers_replications=tillers_replications,
            update_parameters_all_models=RERmax_vegetative_stages_example,
            SENESCWHEAT_TIMESTEP=senescwheat_timestep,
            LIGHT_TIMESTEP=light_timestep,
        )
        positions_greed['wheat'+(i+1)] = ((i+cell_size)%8,(i+cell_size)//8,0) #update for genericity

    lighting = Light_wrapper(
        lightmodel="caribu", 
        out_folder=out_folder, 
        sky=sky,
        planter=planter, 
        indexer=index_log,
        writegeo=write_geo
    )

    current_time_of_the_system = time.time()

    if start_wheat is not None: 
        wheat_dict["wheat1"].start_time=wheat_dict["wheat1"].meteo[wheat_dict["wheat1"].meteo['Date']==start_wheat].index[0]

    for t in range(wheat_dict["wheat1"].start_time, wheat_dict["wheat1"].start_time+simulation_length, wheat_dict["wheat1"].SENESCWHEAT_TIMESTEP):
        if (t % light_timestep == 0) and (wheat_dict["wheat1"].PARi_next_hours(t) > 0):
            
            for i in range (0,n_wheat):

                wheat_input, stems = wheat_dict['wheat'+(i+1)].light_inputs(planter)
                lighting.run(scenes=[wheat_input], day=wheat_dict['wheat'+(i+1)].doy(t), hour=wheat_dict['wheat'+(i+1)].hour(t), parunit="micromol.m-2.s-1", stems=stems)
                wheat_dict['wheat'+(i+1)].light_results(energy=wheat_dict['wheat'+(i+1)].energy(t), lighting=lighting)

        wheat_dict['wheat'+(i+1)].run(t)

    execution_time = int(time.time() - current_time_of_the_system)
    print("\n" "Simulation run in {}".format(str(datetime.timedelta(seconds=execution_time))))


    wheat_dict['wheat'+(i+1)].end(run_postprocessing=run_postprocessing,run_graphs=run_graphs)


if __name__ == "__main__":
    in_folder = "inputs_fspmwheat"
    out_folder = "outputs/cnwheat_default"
    start_wheat = None
    simulation_length = 2500
    write_geo = True
    run_postprocessing=False
    run_graphs=False

    simulation(in_folder, out_folder,start_wheat, simulation_length, write_geo=write_geo,run_postprocessing=run_postprocessing,run_graphs=run_graphs)
