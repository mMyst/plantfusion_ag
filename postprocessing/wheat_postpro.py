from plantfusion.wheat_wrapper import Wheat_wrapper
from plantfusion.indexer import Indexer
from plantfusion.planter import Planter

def w_postpro(out_folder, run_postprocessing=True, run_graphs=True):
    in_folder = r'C:\Users\agrumel\Code\Python_Ecophy\plantfusion_ag\inputs_fspmwheat'
    plants_name = "wheat"
    index_log = Indexer(global_order=[plants_name], wheat_names=[plants_name])
    N_fertilizations = {2016: 357143, 2520: 1000000}
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

    planter = Planter(generation_type="default", indexer=index_log, inter_rows=0.15, plant_density=plant_density)

    wheat = Wheat_wrapper(
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

    wheat.end(simu_ran=False, run_postprocessing=run_postprocessing, run_graphs=run_graphs)

if __name__ == "__main__":
    out_folder = r'C:\Users\agrumel\Code\Python_Ecophy\plantfusion_ag\outputs\cnwheat_soil3ds\bound\0.2m'
    w_postpro(out_folder, run_postprocessing=True, run_graphs=True)
