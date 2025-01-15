
from postprocessing.wheat_postpro import w_postpro 


# Define the range of N values to test
N_values = [0,30,50,100,160,250,350,450]

# Define a function to run the simulation for a single N value
    
def run_postprocessings(N):
    w_postpro(out_folder= "outputs/cnwheat_soil3ds_" + str(N)+'N',run_postprocessing=False,run_graphs=True)

if __name__ == '__main__':
    # Loop through each N value and run the postprocessing
    for N in N_values:
        run_postprocessings(N)
