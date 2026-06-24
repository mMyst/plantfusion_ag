from plantfusion.new_wheat_wrapper import Wheat_wrapper
from plantfusion.light_wrapper import Light_wrapper
from plantfusion.indexer import Indexer
from plantfusion.planter import Planter

import time
import datetime
import os


def simulation(in_folder, out_folder,start_wheat=None, simulation_length=4000, write_geo=False, run_postprocessing=False, run_graphs=False):
    try:
        # Create target Directory
        os.mkdir(os.path.normpath(out_folder))
        print("Directory ", os.path.normpath(out_folder), " Created ")
    except FileExistsError:
        print("Directory ", os.path.normpath(out_folder), " already exists")

    plants_name = "wheat"
    index_log = Indexer(global_order=[plants_name], wheat_names=[plants_name])


    N_fertilizations = {2949: 357143, 4029: 1000000} #22 février et 1er avril
    tillers_replications = {"T1": 0.5, "T2": 0.5, "T3": 0.5, "T4": 0.5}
    plant_density = {1: 250}
    sky = "turtle46" #[4, 5, "soc"]
    RERmax_vegetative_stages_example = {
        "elongwheat": {
            "RERmax": {5: 3.35e-06, 6: 2.1e-06, 7: 2.0e-06, 8: 1.83e-06, 9: 1.8e-06, 10: 1.65e-06, 11: 1.56e-06}
        }
    }
    senescwheat_timestep = 1
    light_timestep = 4
    planter = Planter(generation_type="default", indexer=index_log, inter_rows=0.15, plant_density=plant_density)

    # RATP parameters
    dv = 0.05
    voxels_size = [dv, dv, dv]

    wheat_ratp = Wheat_wrapper(
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

    lighting_ratp = Light_wrapper(
        lightmodel="ratp", 
        out_folder=out_folder, 
        sky=sky,
        planter=planter, 
        indexer=index_log,
        voxels_size=voxels_size,
        angle_distrib_algo="compute global",
        writegeo=write_geo
    )

    light_data = {"PARa": [], "t": []}

    ### SIMULATION ###


    try:
        current_time_of_the_system = time.time()
        for t in range(wheat_ratp.start_time, simulation_length, wheat_ratp.SENESCWHEAT_TIMESTEP):
            if (t % light_timestep == 0) and (wheat_ratp.PARi_next_hours(t) > 0) and bool(wheat_ratp.g.property('geometry')):
                wheat_input, stems = wheat_ratp.light_inputs(planter)
                
                start = time.time()

                lighting_ratp.run(scenes=[wheat_input], day=wheat_ratp.doy(t), hour=wheat_ratp.hour(t), parunit="micromol.m-2.s-1", stems=stems)
                ratp_time = time.time() - start
                
                wheat_ratp.light_results(energy=wheat_ratp.energy (t), lighting=lighting_ratp)

                print("Lighting running time | RATP: ",ratp_time)
            
            print("Step "+ str(t-wheat_ratp.start_time)+"/"+str(simulation_length))
            wheat_ratp.run(t)

        execution_time = int(time.time() - current_time_of_the_system)
        print("\n" "Simulation run in {}".format(str(datetime.timedelta(seconds=execution_time))))

    finally:
        wheat_ratp.end(run_postprocessing=run_postprocessing, run_graphs=run_graphs)


if __name__ == "__main__":
    in_folder = "inputs_fspmwheat"
    out_folder = "outputs/cnwheat_ratp_emergence"
    start_wheat = None
    simulation_length = 4000
    write_geo = True
    run_postprocessing=True
    run_graphs=True

    simulation(in_folder, out_folder,start_wheat, simulation_length, write_geo=write_geo,run_postprocessing=run_postprocessing,run_graphs=run_graphs)
