import concurrent.futures
import os
from simulations.wheat_soil3ds_rangeN import simulation

from postprocessing.wheat_postpro import w_postpro 

# Define the range of N values to test
N_values = [0,30,50,100,160,250,350,450]

in_folder = "inputs_fspmwheat"

start_wheat = None
simulation_length = 2500
write_geo = True

# Define a function to run the simulation for a single N value

def run_simulation(N):
    simulation(in_folder, 
               out_folder="outputs/cnwheat_soil3ds" + '_'+ str(N)+'N', 
               idusm=N, 
               start_wheat=start_wheat, 
               simulation_length=simulation_length,
               write_geo=write_geo)
    
def run_postprocessings(N):
    w_postpro(out_folder= r'C:\Users\agrumel\Documents\Données\Sorties CNWheat\soil3DS_init60N_v2\cnwheat_soil3ds_' + str(N)+'N',run_postprocessing=True,run_graphs=True)

if __name__ == '__main__':
    # Get the number of available cores
    num_cores = os.cpu_count()

    # Use all cores except one
    num_workers = num_cores - 1

    # Run the simulations for each N value simultaneously on all cores except one
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_workers) as executor:
        executor.map(run_postprocessings, N_values)
