import concurrent.futures
import os
from simulations.cn_wheat.wheat_soil3ds_rangeN import simulation
from postprocessing.wheat_postpro import w_postpro

# Define the range of N values to test
N_values = [0, 30, 50, 100, 160, 250, 350, 450]

in_folder = "inputs_fspmwheat"
start_wheat = None
simulation_length = 2500
write_geo = True
sim_rtype = "profile"
sim_min_depth = 0.20  # m
usmfolder = 'lowN' #mixN or lowN


#sim_rtype = "homogeneous"
#sim_min_depth = "max_"

# Define a function to run the simulation for a single N value
def run_simulation(N):
    try:
        sim_out_folder = f"outputs/cnwheat_soil3ds/{sim_rtype}/{usmfolder}/{sim_min_depth}m_{N}N"
        print(f"Running simulation for N={N} with output folder: {sim_out_folder}")

        simulation(
            in_folder,
            out_folder=sim_out_folder,
            rootdistribtype=sim_rtype,  # "bound" or "homogeneous"
            min_depth=sim_min_depth,  # m
            idusm=N,
            start_wheat=start_wheat,
            simulation_length=simulation_length,
            write_geo=write_geo,
            geostep=100,
            run_postprocessing=True,
            run_graphs=True,
            usmfolder=usmfolder
        )
        print(f"Simulation for N={N} completed successfully.")
    except Exception as e:
        print(f"Error running simulation for N={N}: {e}")

def run_postprocessings(N):
    try:
        out_folder = f"C:\\Users\\agrumel\\Documents\\Données\\Sorties CNWheat\\soil3DS_roots\\cnwheat_soil3ds_{N}N"
        print(f"Running postprocessing for N={N} with output folder: {out_folder}")

        w_postpro(
            out_folder=out_folder,
            run_postprocessing=True,
            run_graphs=True
        )
        print(f"Postprocessing for N={N} completed successfully.")
    except Exception as e:
        print(f"Error running postprocessing for N={N}: {e}")

if __name__ == '__main__':
    # Get the number of available cores
    num_cores = os.cpu_count()

    # Use all cores except one
    num_workers = num_cores - 5
    print(f"Using {num_workers} workers out of {num_cores} available cores.")

    # Run the simulations for each N value simultaneously on all cores except one
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_workers) as executor:
        executor.map(run_simulation, N_values)

