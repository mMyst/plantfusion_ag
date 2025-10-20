"""

    contains Wheat_wrapper

"""

import os
import numpy
import pandas
import math
import statsmodels.api as sm

from alinea.adel.adel_dynamic import AdelDyn
from alinea.adel.echap_leaf import echap_leaves

from lightvegemanager.stems import extract_stems_from_MTG

from elongwheat import parameters as elongwheat_parameters

from fspmwheat import caribu_facade
from fspmwheat import cnwheat_facade
from fspmwheat import elongwheat_facade
from fspmwheat import farquharwheat_facade
from fspmwheat import growthwheat_facade
from fspmwheat import senescwheat_facade
from fspmwheat import fspmwheat_facade

from cnwheat import (
    simulation as cnwheat_simulation,
    model as cnwheat_model,
)

from soil3ds.IOxls import read_plant_param

from plantfusion.utils import create_child_folder, save_df_to_csv
from plantfusion.planter import Planter
from plantfusion.indexer import Indexer

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick


class Wheat_wrapper(object):
    """Wrapper for WheatFspm model

    This FSPM is composed with many submodels.

    It creates and stores an instance of
    * AdelDyn
    * openalea.MTG
    * fspmwheat_facade
    * caribu_facade
    * cnwheat_facade
    * elongwheat_facade
    * farquharwheat_facade
    * growthwheat_facade
    * senescwheat_facade

    It also stores all file path names. There is the possibility to restart simulation from previous data.

    Parameters
    ----------
    name : str, optional
        name of the fspm instance, by default "wheat"
    in_folder : str, optional
        input folder path, by default ""
    out_folder : str, optional
        output folder path if writing outputs is activated, by default None
    N_fertilizations : dict, optional
        when to fertilize and how many quantity, {timestep : mol N}, by default {}
    tillers_replications : dict, optional
        triggers tillers replication {tiller tag : length}, by default {}
    planter : Planter, optional
        Object containing plant positions and/or number of plants and/or soil domain, by default None
    indexer : Indexer, optional
        indexer for listing FSPM in the simulation, by default Indexer()
    run_from_outputs : bool, optional
        if you want to start simulation from previous results, by default False
    external_soil_model : bool, optional
        deactivate internal soil process in cn-wheat, by default False
    nitrates_uptake_forced : bool, optional
        forced nitrates uptake data, only external_soil_model is True, by default False
    initialize_nitrates_uptake : float, optional
        N uptake before the simulation start, by default 0.25
    update_parameters_all_models : dict, optional
        forced model parameters, by default None
    stored_times : str or list or None, optional
        stored all or specific simulation timestep, by default None
    option_static : bool, optional
        _description_, by default False
    LIGHT_TIMESTEP : int, optional
        number of iteration between two lighting step, by default 4
    SENESCWHEAT_TIMESTEP : int, optional
        number of iteration between two senescence step, by default 1
    FARQUHARWHEAT_TIMESTEP : int, optional
        number of iteration between two farquhar step, by default 1
    ELONGWHEAT_TIMESTEP : int, optional
        number of iteration between two elongwheat step, by default 1
    GROWTHWHEAT_TIMESTEP : int, optional
        number of iteration between two growthwheat step, by default 1
    CNWHEAT_TIMESTEP : int, optional
        number of iteration between two cnwheat step, by default 1
    AXES_INITIAL_STATE_FILENAME : str, optional
        initialization file name for axes, by default "axes_initial_state.csv"
    ORGANS_INITIAL_STATE_FILENAME : str, optional
        initialization file name for organs, by default "organs_initial_state.csv"
    HIDDENZONES_INITIAL_STATE_FILENAME : str, optional
        initialization file name for hiddenzones, by default "hiddenzones_initial_state.csv"
    ELEMENTS_INITIAL_STATE_FILENAME : str, optional
        initialization file name for elements, by default "elements_initial_state.csv"
    SOILS_INITIAL_STATE_FILENAME : str, optional
        initialization file name for soil, by default "soils_initial_state.csv"
    METEO_FILENAME : str, optional
        meteo file name, by default "meteo_Ljutovac2002.csv"
    NITRATES_UPTAKE_FORCINGS_FILENAME : str, optional
        forced nitrates uptake file name, by default "nitrates_uptake_forcings.csv"
    SOIL_PARAMETERS_FILENAME : str, optional
        WheatFspm soil parameters for soil3ds, by default ""
    SOIL_PARAMETERS_SHEETNAME : str, optional
        sheet in SOIL_PARAMETERS_FILENAME corresponding current plant specy parameters, by default "cnwheat"
    AXES_OUTPUTS_FILENAME : str, optional
        brut outputs file name for axes, by default "axes_outputs.csv"
    ORGANS_OUTPUTS_FILENAME : str, optional
        brut outputs file name for organs, by default "organs_outputs.csv"
    HIDDENZONES_OUTPUTS_FILENAME : str, optional
        brut outputs file name for hiddenzones, by default "hiddenzones_outputs.csv"
    ELEMENTS_OUTPUTS_FILENAME : str, optional
        brut outputs file name for elements, by default "elements_outputs.csv"
    SOILS_OUTPUTS_FILENAME : str, optional
        brut outputs file name for soil, by default "soils_outputs.csv"
    AXES_POSTPROCESSING_FILENAME : str, optional
        postprocessing outputs file name for axes, by default "axes_postprocessing.csv"
    ORGANS_POSTPROCESSING_FILENAME : str, optional
        postprocessing outputs file name for organs, by default "organs_postprocessing.csv"
    HIDDENZONES_POSTPROCESSING_FILENAME : str, optional
        postprocessing outputs file name for hiddenzones, by default "hiddenzones_postprocessing.csv"
    ELEMENTS_POSTPROCESSING_FILENAME : str, optional
        postprocessing outputs file name for elements, by default "elements_postprocessing.csv"
    SOILS_POSTPROCESSING_FILENAME : str, optional
        postprocessing outputs file name for soil, by default "soils_postprocessing.csv"

    """    

    # constants in class
    HOUR_TO_SECOND_CONVERSION_FACTOR = 3600
    AXES_INDEX_COLUMNS = ["t", "plant", "axis"]
    ELEMENTS_INDEX_COLUMNS = ["t", "plant", "axis", "metamer", "organ", "element"]
    HIDDENZONES_INDEX_COLUMNS = ["t", "plant", "axis", "metamer"]
    ORGANS_INDEX_COLUMNS = ["t", "plant", "axis", "organ"]
    SOILS_INDEX_COLUMNS = ["t", "plant", "axis"]

    def __init__(
        self,
        name="wheat",
        in_folder="",
        out_folder=None,
        N_fertilizations={},
        tillers_replications={},
        planter=Planter(),
        indexer=Indexer(),
        run_from_outputs=False,
        external_soil_model=False,
        nitrates_uptake_forced=False,
        initialize_nitrates_uptake=0.25,
        update_parameters_all_models=None,
        stored_times=None,
        option_static=False,
        LIGHT_TIMESTEP=4,
        SENESCWHEAT_TIMESTEP=1,
        FARQUHARWHEAT_TIMESTEP=1,
        ELONGWHEAT_TIMESTEP=1,
        GROWTHWHEAT_TIMESTEP=1,
        CNWHEAT_TIMESTEP=1,
        AXES_INITIAL_STATE_FILENAME="axes_initial_state.csv",
        ORGANS_INITIAL_STATE_FILENAME="organs_initial_state.csv",
        HIDDENZONES_INITIAL_STATE_FILENAME="hiddenzones_initial_state.csv",
        ELEMENTS_INITIAL_STATE_FILENAME="elements_initial_state.csv",
        SOILS_INITIAL_STATE_FILENAME="soils_initial_state.csv",
        METEO_FILENAME="meteo_Ljutovac2002.csv",
        NITRATES_UPTAKE_FORCINGS_FILENAME="nitrates_uptake_forcings.csv",
        SOIL_PARAMETERS_FILENAME="",
        SOIL_PARAMETERS_SHEETNAME="cnwheat",
        AXES_OUTPUTS_FILENAME="axes_outputs.csv",
        ORGANS_OUTPUTS_FILENAME="organs_outputs.csv",
        HIDDENZONES_OUTPUTS_FILENAME="hiddenzones_outputs.csv",
        ELEMENTS_OUTPUTS_FILENAME="elements_outputs.csv",
        SOILS_OUTPUTS_FILENAME="soils_outputs.csv",
        AXES_POSTPROCESSING_FILENAME="axes_postprocessing.csv",
        ORGANS_POSTPROCESSING_FILENAME="organs_postprocessing.csv",
        HIDDENZONES_POSTPROCESSING_FILENAME="hiddenzones_postprocessing.csv",
        ELEMENTS_POSTPROCESSING_FILENAME="elements_postprocessing.csv",
        SOILS_POSTPROCESSING_FILENAME="soils_postprocessing.csv",
        rootdistribtype = "homogeneous"
    ) -> None:
        """Constructor, 
    
        """ 

        self.N_fertilizations = N_fertilizations
        self.tillers_replications = tillers_replications

        self.external_soil_model = external_soil_model
        self.nitrates_uptake_forced = nitrates_uptake_forced
        self.option_static = option_static

        self.last_year_doy = 0

        self.name = name
        self.indexer = indexer
        self.global_index = indexer.global_order.index(name)
        self.wheat_index = indexer.wheat_names.index(name)

        self.nb_plants = planter.number_of_plants[self.global_index]
        if name in planter.plant_density:
            self.plant_density = {1: planter.plant_density[name]}
        else:
            self.plant_density = planter.plant_density
        self.generation_type = planter.generation_type

        self.LIGHT_TIMESTEP = LIGHT_TIMESTEP
        self.SENESCWHEAT_TIMESTEP = SENESCWHEAT_TIMESTEP
        self.FARQUHARWHEAT_TIMESTEP = FARQUHARWHEAT_TIMESTEP
        self.ELONGWHEAT_TIMESTEP = ELONGWHEAT_TIMESTEP
        self.GROWTHWHEAT_TIMESTEP = GROWTHWHEAT_TIMESTEP
        self.CNWHEAT_TIMESTEP = CNWHEAT_TIMESTEP

        self.AXES_OUTPUTS_FILENAME = AXES_OUTPUTS_FILENAME
        self.ORGANS_OUTPUTS_FILENAME = ORGANS_OUTPUTS_FILENAME
        self.HIDDENZONES_OUTPUTS_FILENAME = HIDDENZONES_OUTPUTS_FILENAME
        self.ELEMENTS_OUTPUTS_FILENAME = ELEMENTS_OUTPUTS_FILENAME
        self.SOILS_OUTPUTS_FILENAME = SOILS_OUTPUTS_FILENAME

        in_folder = os.path.normpath(in_folder)

        self.AXES_POSTPROCESSING_FILENAME = AXES_POSTPROCESSING_FILENAME
        self.ORGANS_POSTPROCESSING_FILENAME = ORGANS_POSTPROCESSING_FILENAME
        self.HIDDENZONES_POSTPROCESSING_FILENAME = HIDDENZONES_POSTPROCESSING_FILENAME
        self.ELEMENTS_POSTPROCESSING_FILENAME = ELEMENTS_POSTPROCESSING_FILENAME
        self.SOILS_POSTPROCESSING_FILENAME = SOILS_POSTPROCESSING_FILENAME

        if out_folder is not None:
            self.out_folder = os.path.join(os.path.normpath(out_folder), name)
            try:
                os.mkdir(os.path.normpath(self.out_folder))
                print("Directory ", self.out_folder, " Created ")
            except FileExistsError:
                pass

            create_child_folder(self.out_folder, "brut")
            create_child_folder(self.out_folder, "postprocessing")
            create_child_folder(self.out_folder, "graphs")

        if external_soil_model:
            SOILS_INITIAL_STATE_FILENAME = None
            self.rootdistribtype = rootdistribtype
            
            if nitrates_uptake_forced:
                nitrates_uptake_data_filepath = os.path.join(in_folder, NITRATES_UPTAKE_FORCINGS_FILENAME)
                nitrates_uptake_data_df = pandas.read_csv(nitrates_uptake_data_filepath)
                self.nitrates_uptake_data_grouped = nitrates_uptake_data_df.groupby(
                    cnwheat_simulation.Simulation.ORGANS_T_INDEXES
                )

            else:
                wheat_paramp = read_plant_param(os.path.normpath(SOIL_PARAMETERS_FILENAME), SOIL_PARAMETERS_SHEETNAME)
                self.soil_parameters = [wheat_paramp] * self.nb_plants

        self.inputs_dataframes = {}
        new_start_time = -1
        if run_from_outputs:
            previous_outputs_dataframes = {}

            for initial_state_filename, outputs_filename, index_columns in (
                (AXES_INITIAL_STATE_FILENAME, AXES_OUTPUTS_FILENAME, self.AXES_INDEX_COLUMNS),
                (ORGANS_INITIAL_STATE_FILENAME, ORGANS_OUTPUTS_FILENAME, self.ORGANS_INDEX_COLUMNS),
                (HIDDENZONES_INITIAL_STATE_FILENAME, HIDDENZONES_OUTPUTS_FILENAME, self.HIDDENZONES_INDEX_COLUMNS),
                (ELEMENTS_INITIAL_STATE_FILENAME, ELEMENTS_OUTPUTS_FILENAME, self.ELEMENTS_INDEX_COLUMNS),
                (SOILS_INITIAL_STATE_FILENAME, SOILS_OUTPUTS_FILENAME, self.SOILS_INDEX_COLUMNS),
            ):
                if not initial_state_filename:
                    continue
                previous_outputs_dataframe = pandas.read_csv(os.path.join(out_folder, outputs_filename))
                # Convert NaN to None
                previous_outputs_dataframes[outputs_filename] = previous_outputs_dataframe.replace({numpy.nan: None})

                assert "t" in previous_outputs_dataframes[outputs_filename].columns

                last_t_step = max(previous_outputs_dataframes[outputs_filename]["t"])
                new_start_time = last_t_step + 1

                if initial_state_filename == ELEMENTS_INITIAL_STATE_FILENAME:
                    elements_previous_outputs = previous_outputs_dataframes[outputs_filename]
                    new_initial_state = elements_previous_outputs[~elements_previous_outputs.is_over.isnull()]
                else:
                    new_initial_state = previous_outputs_dataframes[outputs_filename]
                idx = (
                    new_initial_state.groupby([col for col in index_columns if col != "t"])["t"].transform(max)
                    == new_initial_state["t"]
                )
                self.inputs_dataframes[initial_state_filename] = new_initial_state[idx].drop(["t"], axis=1)

            # Make sure boolean columns have either type bool or float
            bool_columns = [
                "is_over",
                "is_growing",
                "leaf_is_emerged",
                "internode_is_visible",
                "leaf_is_growing",
                "internode_is_growing",
                "leaf_is_remobilizing",
                "internode_is_remobilizing",
            ]
            for df in [
                self.inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME],
                self.inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME],
            ]:
                for cln in bool_columns:
                    if cln in df.keys():
                        df[cln].replace(to_replace="False", value=0.0, inplace=True)
                        df[cln].replace(to_replace="True", value=1.0, inplace=True)
                        df[cln] = pandas.to_numeric(df[cln])
        else:
            for inputs_filename in (
                AXES_INITIAL_STATE_FILENAME,
                ORGANS_INITIAL_STATE_FILENAME,
                HIDDENZONES_INITIAL_STATE_FILENAME,
                ELEMENTS_INITIAL_STATE_FILENAME,
                SOILS_INITIAL_STATE_FILENAME,
            ):
                if not inputs_filename:
                    continue
                inputs_dataframe = pandas.read_csv(os.path.join(in_folder, inputs_filename))
                self.inputs_dataframes[inputs_filename] = inputs_dataframe.replace({numpy.nan: None})

        # Start time of the simulation
        self.start_time = max(0, new_start_time)

        # Name of the CSV files which contains the meteo data
        self.meteo = pandas.read_csv(os.path.join(in_folder, METEO_FILENAME), index_col="t")

        # -- OUTPUTS CONFIGURATION --

        # Save the outputs with a full scan of the MTG at each time step (or at selected time steps)
        UPDATE_SHARED_DF = False
        if stored_times is None:
            stored_times = "all"
        if not (stored_times == "all" or isinstance(stored_times, list)):
            print("stored_times should be either 'all', a list or an empty list.")
            raise
        self.stored_times = stored_times

        # create empty dataframes to shared data between the models
        self.shared_axes_inputs_outputs_df = pandas.DataFrame()
        self.shared_organs_inputs_outputs_df = pandas.DataFrame()
        self.shared_hiddenzones_inputs_outputs_df = pandas.DataFrame()
        self.shared_elements_inputs_outputs_df = pandas.DataFrame()
        self.shared_soils_inputs_outputs_df = pandas.DataFrame()

        # define lists of dataframes to store the inputs and the outputs of the models at each step.
        self.axes_all_data_list = []
        self.organs_all_data_list = []
        self.hiddenzones_all_data_list = []
        self.elements_all_data_list = []
        self.soils_all_data_list = []

        self.all_simulation_steps = []

        # -- ADEL and MTG CONFIGURATION --

        # read adelwheat inputs at t0
        self.adel_wheat = AdelDyn(seed=1, scene_unit="m", leaves=echap_leaves(xy_model="Soissons_byleafclass"))
        self.g = self.adel_wheat.load(dir=in_folder)

        # ---------------------------------------------
        # ----- CONFIGURATION OF THE wrapperS -------
        # ---------------------------------------------

        # -- ELONGWHEAT (created first because it is the only wrapper to add new metamers) --
        # Initial states
        elongwheat_hiddenzones_initial_state = self.inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME][
            elongwheat_facade.converter.HIDDENZONE_TOPOLOGY_COLUMNS
            + [
                i
                for i in elongwheat_facade.simulation.HIDDENZONE_INPUTS
                if i in self.inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME].columns
            ]
        ].copy()
        elongwheat_elements_initial_state = self.inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME][
            elongwheat_facade.converter.ELEMENT_TOPOLOGY_COLUMNS
            + [
                i
                for i in elongwheat_facade.simulation.ELEMENT_INPUTS
                if i in self.inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME].columns
            ]
        ].copy()
        elongwheat_axes_initial_state = self.inputs_dataframes[AXES_INITIAL_STATE_FILENAME][
            elongwheat_facade.converter.AXIS_TOPOLOGY_COLUMNS
            + [
                i
                for i in elongwheat_facade.simulation.AXIS_INPUTS
                if i in self.inputs_dataframes[AXES_INITIAL_STATE_FILENAME].columns
            ]
        ].copy()

        phytoT_df = pandas.read_csv(os.path.join(in_folder, "phytoT.csv"))

        # update parameters if specified
        if update_parameters_all_models and "elongwheat" in update_parameters_all_models:
            update_parameters_elongwheat = update_parameters_all_models["elongwheat"]
        else:
            update_parameters_elongwheat = None

        # wrapper initialisation
        self.elongwheat_facade_ = elongwheat_facade.ElongWheatFacade(
            self.g,
            ELONGWHEAT_TIMESTEP * self.HOUR_TO_SECOND_CONVERSION_FACTOR,
            elongwheat_axes_initial_state,
            elongwheat_hiddenzones_initial_state,
            elongwheat_elements_initial_state,
            self.shared_axes_inputs_outputs_df,
            self.shared_hiddenzones_inputs_outputs_df,
            self.shared_elements_inputs_outputs_df,
            self.adel_wheat,
            phytoT_df,
            update_parameters_elongwheat,
            update_shared_df=UPDATE_SHARED_DF,
        )

        # -- CARIBU --
        self.caribu_facade_ = caribu_facade.CaribuFacade(
            self.g, self.shared_elements_inputs_outputs_df, self.adel_wheat, update_shared_df=UPDATE_SHARED_DF
        )

        # -- SENESCWHEAT --
        # Initial states
        senescwheat_roots_initial_state = (
            self.inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME]
            .loc[self.inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME]["organ"] == "roots"][
                senescwheat_facade.converter.ROOTS_TOPOLOGY_COLUMNS
                + [
                    i
                    for i in senescwheat_facade.converter.SENESCWHEAT_ROOTS_INPUTS
                    if i in self.inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME].columns
                ]
            ]
            .copy()
        )

        senescwheat_elements_initial_state = self.inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME][
            senescwheat_facade.converter.ELEMENTS_TOPOLOGY_COLUMNS
            + [
                i
                for i in senescwheat_facade.converter.SENESCWHEAT_ELEMENTS_INPUTS
                if i in self.inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME].columns
            ]
        ].copy()

        senescwheat_axes_initial_state = self.inputs_dataframes[AXES_INITIAL_STATE_FILENAME][
            senescwheat_facade.converter.AXES_TOPOLOGY_COLUMNS
            + [
                i
                for i in senescwheat_facade.converter.SENESCWHEAT_AXES_INPUTS
                if i in self.inputs_dataframes[AXES_INITIAL_STATE_FILENAME].columns
            ]
        ].copy()

        # update parameters if specified
        if update_parameters_all_models and "senescwheat" in update_parameters_all_models:
            update_parameters_senescwheat = update_parameters_all_models["senescwheat"]
        else:
            update_parameters_senescwheat = None

        # wrapper initialisation
        self.senescwheat_facade_ = senescwheat_facade.SenescWheatFacade(
            self.g,
            SENESCWHEAT_TIMESTEP * self.HOUR_TO_SECOND_CONVERSION_FACTOR,
            senescwheat_roots_initial_state,
            senescwheat_axes_initial_state,
            senescwheat_elements_initial_state,
            self.shared_organs_inputs_outputs_df,
            self.shared_axes_inputs_outputs_df,
            self.shared_elements_inputs_outputs_df,
            update_parameters_senescwheat,
            update_shared_df=UPDATE_SHARED_DF,
        )

        # -- FARQUHARWHEAT --
        # Initial states
        farquharwheat_elements_initial_state = self.inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME][
            farquharwheat_facade.converter.ELEMENT_TOPOLOGY_COLUMNS
            + [
                i
                for i in farquharwheat_facade.converter.FARQUHARWHEAT_ELEMENTS_INPUTS
                if i in self.inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME].columns
            ]
        ].copy()

        farquharwheat_axes_initial_state = self.inputs_dataframes[AXES_INITIAL_STATE_FILENAME][
            farquharwheat_facade.converter.AXIS_TOPOLOGY_COLUMNS
            + [
                i
                for i in farquharwheat_facade.converter.FARQUHARWHEAT_AXES_INPUTS
                if i in self.inputs_dataframes[AXES_INITIAL_STATE_FILENAME].columns
            ]
        ].copy()

        # Use the initial version of the photosynthesis sub-model (as in Barillot et al. 2016, and in Gauthier et al. 2020)
        update_parameters_farquharwheat = {"SurfacicProteins": False, "NSC_Retroinhibition": False}

        # wrapper initialisation
        self.farquharwheat_facade_ = farquharwheat_facade.FarquharWheatFacade(
            self.g,
            farquharwheat_elements_initial_state,
            farquharwheat_axes_initial_state,
            self.shared_elements_inputs_outputs_df,
            update_parameters_farquharwheat,
            update_shared_df=UPDATE_SHARED_DF,
        )

        # -- GROWTHWHEAT --
        # Initial states
        growthwheat_hiddenzones_initial_state = self.inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME][
            growthwheat_facade.converter.HIDDENZONE_TOPOLOGY_COLUMNS
            + [
                i
                for i in growthwheat_facade.simulation.HIDDENZONE_INPUTS
                if i in self.inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME].columns
            ]
        ].copy()

        growthwheat_elements_initial_state = self.inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME][
            growthwheat_facade.converter.ELEMENT_TOPOLOGY_COLUMNS
            + [
                i
                for i in growthwheat_facade.simulation.ELEMENT_INPUTS
                if i in self.inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME].columns
            ]
        ].copy()

        growthwheat_root_initial_state = (
            self.inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME]
            .loc[self.inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME]["organ"] == "roots"][
                growthwheat_facade.converter.ROOT_TOPOLOGY_COLUMNS
                + [
                    i
                    for i in growthwheat_facade.simulation.ROOT_INPUTS
                    if i in self.inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME].columns
                ]
            ]
            .copy()
        )

        growthwheat_axes_initial_state = self.inputs_dataframes[AXES_INITIAL_STATE_FILENAME][
            growthwheat_facade.converter.AXIS_TOPOLOGY_COLUMNS
            + [
                i
                for i in growthwheat_facade.simulation.AXIS_INPUTS
                if i in self.inputs_dataframes[AXES_INITIAL_STATE_FILENAME].columns
            ]
        ].copy()

        # update parameters if specified
        if update_parameters_all_models and "growthwheat" in update_parameters_all_models:
            update_parameters_growthwheat = update_parameters_all_models["growthwheat"]
        else:
            update_parameters_growthwheat = None

        # wrapper initialisation
        self.growthwheat_facade_ = growthwheat_facade.GrowthWheatFacade(
            self.g,
            GROWTHWHEAT_TIMESTEP * self.HOUR_TO_SECOND_CONVERSION_FACTOR,
            growthwheat_hiddenzones_initial_state,
            growthwheat_elements_initial_state,
            growthwheat_root_initial_state,
            growthwheat_axes_initial_state,
            self.shared_organs_inputs_outputs_df,
            self.shared_hiddenzones_inputs_outputs_df,
            self.shared_elements_inputs_outputs_df,
            self.shared_axes_inputs_outputs_df,
            update_parameters_growthwheat,
            update_shared_df=UPDATE_SHARED_DF,
        )

        # -- CNWHEAT --
        # Initial states
        cnwheat_organs_initial_state = self.inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME][
            [
                i
                for i in cnwheat_facade.cnwheat_converter.ORGANS_VARIABLES
                if i in self.inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME].columns
            ]
        ].copy()

        cnwheat_hiddenzones_initial_state = self.inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME][
            [
                i
                for i in cnwheat_facade.cnwheat_converter.HIDDENZONE_VARIABLES
                if i in self.inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME].columns
            ]
        ].copy()

        cnwheat_elements_initial_state = self.inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME][
            [
                i
                for i in cnwheat_facade.cnwheat_converter.ELEMENTS_VARIABLES
                if i in self.inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME].columns
            ]
        ].copy()

        if not external_soil_model:
            cnwheat_soils_initial_state = self.inputs_dataframes[SOILS_INITIAL_STATE_FILENAME][
                [
                    i
                    for i in cnwheat_facade.cnwheat_converter.SOILS_VARIABLES
                    if i in self.inputs_dataframes[SOILS_INITIAL_STATE_FILENAME].columns
                ]
            ].copy()
        else:
            cnwheat_soils_initial_state = None

        # update parameters if specified
        if update_parameters_all_models and "cnwheat" in update_parameters_all_models:
            update_parameters_cnwheat = update_parameters_all_models["cnwheat"]
        else:
            update_parameters_cnwheat = {}

        # wrapper initialisation
        self.cnwheat_facade_ = cnwheat_facade.CNWheatFacade(
            self.g,
            CNWHEAT_TIMESTEP * self.HOUR_TO_SECOND_CONVERSION_FACTOR,
            self.plant_density,
            update_parameters_cnwheat,
            cnwheat_organs_initial_state,
            cnwheat_hiddenzones_initial_state,
            cnwheat_elements_initial_state,
            cnwheat_soils_initial_state,
            self.shared_axes_inputs_outputs_df,
            self.shared_organs_inputs_outputs_df,
            self.shared_hiddenzones_inputs_outputs_df,
            self.shared_elements_inputs_outputs_df,
            self.shared_soils_inputs_outputs_df,
            update_shared_df=UPDATE_SHARED_DF,
            external_soil_model=external_soil_model,
        )

        # Force the senescence and photosynthesis of the population
        if external_soil_model:
            if self.nitrates_uptake_forced:
                t = 0
                self.force_nitrates_uptake(t)
            else:
                self.uptake_nitrate_hour = initialize_nitrates_uptake
                self.update_Nitrates_cnwheat_mtg()

        # -- FSPMWHEAT --
        # wrapper initialisation
        self.fspmwheat_facade_ = fspmwheat_facade.FSPMWheatFacade(self.g)

        # update geometry
        self.adel_wheat.update_geometry(self.g)

    def light_inputs(self, planter):
        """Wheats geometric scene for lighting

        Note
        ----
        It creates a canopy from one wheat geometry with variabilities on each plant

        Parameters
        ----------
        planter : Planter
            provides plant positions and methods for creating the canopy

        Returns
        -------
        plantgl.Scene, list
            
            * wheat canopy from one wheat
            
            * list of tuple precising stems id the plantgl scene. Elements are (specy id, organ id) 
        """        
        if self.generation_type == "default":
            scene_wheat = planter.create_heterogeneous_canopy(
                self.adel_wheat,
                mtg=self.g,
                stem_name="StemElement",
                leaf_name="LeafElement1",
                indice_wheat_instance=self.wheat_index,
            )

        elif self.generation_type == "random":
            scene_wheat = planter.generate_random_wheat(
                self.adel_wheat,
                mtg=self.g,
                indice_wheat_instance=self.wheat_index,
                stem_name="StemElement",
                leaf_name="LeafElement1",
            )

        elif self.generation_type == "row":
            scene_wheat = planter.generate_row_wheat(
                self.adel_wheat, self.g, self.wheat_index, stem_name="StemElement", leaf_name="LeafElement1"
            )
        elif 'forced' in self.generation_type:
            scene_wheat = planter.generate_forced_wheat(
                self.adel_wheat, self.g, self.wheat_index, stem_name="StemElement", leaf_name="LeafElement1"
            )

        else:
            print("can't recognize positions generation type, choose between default, random and row")
            raise

        stems = extract_stems_from_MTG(self.g, self.global_index)

        return scene_wheat, stems

    def light_results(self, energy, lighting, selective_global_index=None):
        """Interprets lighting results

        Parameters
        ----------
        energy : float
            meteo radiation in micromol.m-2.s-1
        lighting : Light_wrapper
            Contains lighting results in a pandas.Dataframe
        selective_global_index : int, optional
            if specy ID in lighting results is different from self.global_index, by default None
        """        
        if selective_global_index is not None:
            saved_global_index = self.global_index
            self.global_index = selective_global_index

        results = lighting.results_organs()
        lightmodel = lighting.lightmodel

        # crée un tableau comme dans caribu_facade de fspm-wheat
        dico_par = {}
        para_dic = {}
        erel_dic = {}

        df_outputs_esp = results[results.VegetationType == self.global_index]
        for s in df_outputs_esp["Organ"]:
            d = df_outputs_esp[df_outputs_esp.Organ == s]

            if lightmodel == "caribu":
                para_dic[s] = d["par Eabs"].values[0] * energy
                erel_dic[s] = d["par Eabs"].values[0]  # lighting ran with energy = 1.

            elif lightmodel == "ratp":
                para_dic[s] = d["PARa"].values[0] * energy
                erel_dic[s] = d["Intercepted"].values[0]

        dico_par["PARa"] = para_dic
        dico_par["Erel"] = erel_dic

        for param in dico_par:
            if param not in self.g.properties():
                self.g.add_property(param)
            # update the self.g
            self.g.property(param).update(dico_par[param])

        if selective_global_index is not None:
            self.global_index = saved_global_index

    def soil_inputs(self, soil, planter, lighting):
        """Compute soil inputs for soil3ds

        Roots length is computed from their mass (g) and a specific root length (m/g).

        Parameters
        ----------
        soil : Soil_wrapper
            for soil3ds dimensions and get plant positions in soil voxels grid
        planter : Planter
            for plant positions
        lighting : Light_wrapper
            for computing light interception capacity per plant

        Returns
        -------
        4-tuple with soil input lists
            
            * (list of float) roots N content per plant [0-1]

            * (numpy.array of 4 dimensions) roots length per plant and soil voxel in m

            * (list of dict) soil parameters per plant (KMAX, VMAX, N content thresholds, ...)

            * (list of float) light capacity interception per plant [0-1]

        """        
        # ls_N
        N_content_roots = self.compute_N_content_roots()
        N_content_roots_per_plant = [N_content_roots] * self.nb_plants

        # ls_roots
        roots_length_per_plant_per_soil_layer = self.compute_roots_length(soil, planter)

        # ls_epsi
        organs_results = lighting.results_organs()
        filtered_data = organs_results[organs_results.VegetationType.isin([self.global_index])]
        plant_leaf_area = numpy.sum(filtered_data["Area"].values) / self.nb_plants
        plants_light_interception = self.compute_plants_light_interception(plant_leaf_area, lighting.soil_energy())

        return (
            N_content_roots_per_plant,
            roots_length_per_plant_per_soil_layer,
            self.soil_parameters,
            plants_light_interception,
        )

    def soil_results(self, uptakeN_per_plant, planter=None, selective_global_index=None):
        """Interprets soil results

        Compute a mean results among all plants.

        Parameters
        ----------
        uptakeN_per_plant : list
            numpy array of 4 dimensions [plant id, nz, nx, ny], N uptake for each plant in each soil voxel
        planter : Planter, optional
            gives index of wheat plants among the result lists, by default None
        selective_global_index : int, optional
            if specy ID in lighting results is different from self.global_index, by default None
        """        
        if selective_global_index is not None:
            saved_global_index = self.global_index
            self.global_index = selective_global_index

        if planter is not None:
            index_in_global_plants = [
                sum(planter.number_of_plants[: self.global_index]),
                sum(planter.number_of_plants[: self.global_index + 1]),
            ]
        else:
            index_in_global_plants = [0, self.nb_plants]

        # considère que tous les blé sont identiques: moyenne
        uptake_nitrate = 0.0
        for i in range(*index_in_global_plants):
            uptake_nitrate += numpy.sum(uptakeN_per_plant[i])
        uptake_nitrate *= 1 / self.nb_plants

        self.uptake_nitrate_hour = self.convert_uptake_nitrate(uptake_nitrate)

        if selective_global_index is not None:
            self.global_index = saved_global_index

    def run(self, t_light):
        """Time step computing of WheatFspm. Independant from other fspm in the simulation

        Parameters
        ----------
        t_light : int
            meteo timestep
        """     

        if not ((t_light % self.LIGHT_TIMESTEP == 0) and (self.PARi_next_hours(t_light) > 0)):
            Erel = self.g.property("Erel")
            PARa_output = {k: v * self.energy(t_light) for k, v in Erel.items()}
            outputs = {}
            outputs.update({"PARa": PARa_output})
            for param in outputs.keys():
                if param not in self.g.properties():
                    self.g.add_property(param)
                # update the MTG
                self.g.property(param).update(outputs[param])

        # suite de la simu
        for t_senescwheat in range(t_light, t_light + self.SENESCWHEAT_TIMESTEP, self.SENESCWHEAT_TIMESTEP):
            # run SenescWheat
            self.senescwheat_facade_.run()

            # Test for dead plant # TODO: adapt in case of multiple plants
            if (
                not self.shared_elements_inputs_outputs_df.empty
                and numpy.nansum(
                    self.shared_elements_inputs_outputs_df.loc[
                        self.shared_elements_inputs_outputs_df["element"].isin(["StemElement", "LeafElement1"]),
                        "green_area",
                    ]
                )
                == 0
            ):
                # append the inputs and outputs at current step to global lists
                self.all_simulation_steps.append(t_senescwheat)
                self.axes_all_data_list.append(self.shared_axes_inputs_outputs_df.copy())
                self.organs_all_data_list.append(self.shared_organs_inputs_outputs_df.copy())
                self.hiddenzones_all_data_list.append(self.shared_hiddenzones_inputs_outputs_df.copy())
                self.elements_all_data_list.append(self.shared_elements_inputs_outputs_df.copy())
                self.soils_all_data_list.append(self.shared_soils_inputs_outputs_df.copy())
                break

            # Run the rest of the model if the plant is alive
            for t_farquharwheat in range(
                t_senescwheat, t_senescwheat + self.SENESCWHEAT_TIMESTEP, self.FARQUHARWHEAT_TIMESTEP
            ):
                # get the meteo of the current step
                Ta, ambient_CO2, RH, Ur = self.meteo.loc[
                    t_farquharwheat, ["air_temperature", "ambient_CO2", "humidity", "Wind"]
                ]

                # run FarquharWheat
                self.farquharwheat_facade_.run(Ta, ambient_CO2, RH, Ur)

                for t_elongwheat in range(
                    t_farquharwheat, t_farquharwheat + self.FARQUHARWHEAT_TIMESTEP, self.ELONGWHEAT_TIMESTEP
                ):
                    # run ElongWheat
                    Tair, Tsoil = self.meteo.loc[t_elongwheat, ["air_temperature", "soil_temperature"]]
                    self.elongwheat_facade_.run(Tair, Tsoil, option_static=self.option_static)

                    # update geometry
                    self.adel_wheat.update_geometry(self.g)

                    for t_growthwheat in range(
                        t_elongwheat, t_elongwheat + self.ELONGWHEAT_TIMESTEP, self.GROWTHWHEAT_TIMESTEP
                    ):
                        # run GrowthWheat
                        self.growthwheat_facade_.run()

                        for t_cnwheat in range(
                            t_growthwheat, t_growthwheat + self.GROWTHWHEAT_TIMESTEP, self.CNWHEAT_TIMESTEP
                        ):
                            print("t cnwheat is {}".format(t_cnwheat))

                            # N fertilization if any
                            if (
                                not self.external_soil_model
                                and self.N_fertilizations is not None
                                and len(self.N_fertilizations) > 0
                            ):
                                if t_cnwheat in self.N_fertilizations.keys():
                                    self.cnwheat_facade_.soils[(1, "MS")].nitrates += self.N_fertilizations[t_cnwheat]

                            if t_cnwheat > 0:
                                if self.external_soil_model:
                                    if self.nitrates_uptake_forced:
                                        self.force_nitrates_uptake(t_cnwheat)
                                    else:
                                        self.update_Nitrates_cnwheat_mtg()

                                # run CNWheat
                                Tair = self.meteo.loc[t_elongwheat, "air_temperature"]
                                Tsoil = self.meteo.loc[t_elongwheat, "soil_temperature"]
                                self.cnwheat_facade_.run(Tair, Tsoil, self.tillers_replications)

                            # append outputs at current step to global lists
                            if (self.stored_times == "all") or (t_cnwheat in self.stored_times):
                                (
                                    axes_outputs,
                                    elements_outputs,
                                    hiddenzones_outputs,
                                    organs_outputs,
                                    soils_outputs,
                                ) = self.fspmwheat_facade_.build_outputs_df_from_MTG()

                                self.all_simulation_steps.append(t_cnwheat)
                                self.axes_all_data_list.append(axes_outputs)
                                self.organs_all_data_list.append(organs_outputs)
                                self.hiddenzones_all_data_list.append(hiddenzones_outputs)
                                self.elements_all_data_list.append(elements_outputs)
                                self.soils_all_data_list.append(soils_outputs)

    def end(self, simu_ran=True, run_postprocessing=False, run_graphs=False):
        """Write output files and close submodel instances

        Parameters
        ----------
        run_postprocessing : bool, optional
            activate postprocessing files writting, by default False
        """        
        PRECISION = 4

        outputs_brut = os.path.join(self.out_folder, "brut")

        if simu_ran==False:#if the simulation has not been ran. Pre existing outputs are imported from outputs_brut
            outputs_df_dict = {}

            for outputs_filename in (self.AXES_OUTPUTS_FILENAME,
                                     self.ORGANS_OUTPUTS_FILENAME,
                                     self.HIDDENZONES_OUTPUTS_FILENAME,
                                     self.ELEMENTS_OUTPUTS_FILENAME,
                                     self.SOILS_OUTPUTS_FILENAME):
                outputs_filepath = os.path.join(outputs_brut, outputs_filename)
                outputs_df = pandas.read_csv(outputs_filepath)
                outputs_file_basename = outputs_filename.split('.')[0]
                outputs_df_dict[outputs_file_basename] = outputs_df

                # Assert states_filepaths were not opened during simulation run meaning that other filenames were saved
                tmp_filename = 'ACTUAL_{}.csv'.format(outputs_file_basename)
                tmp_path = os.path.join(outputs_brut, tmp_filename)
                assert not os.path.isfile(tmp_path), \
                    "File {} was saved because {} was opened during simulation run. Rename it before running postprocessing".format(tmp_filename, outputs_file_basename)

            time_grid = list(outputs_df_dict.values())[0].t
            delta_t = (time_grid.loc[1] - time_grid.loc[0]) * self.HOUR_TO_SECOND_CONVERSION_FACTOR

        else:#if the simulation has been ran, the outputs are written and saved here
            # save des données brutes
            outputs_df_dict = {}
            for outputs_df_list, outputs_filename, index_columns in (
                (self.axes_all_data_list, self.AXES_OUTPUTS_FILENAME, cnwheat_simulation.Simulation.AXES_T_INDEXES),
                (self.organs_all_data_list, self.ORGANS_OUTPUTS_FILENAME, cnwheat_simulation.Simulation.ORGANS_T_INDEXES),
                (
                    self.hiddenzones_all_data_list,
                    self.HIDDENZONES_OUTPUTS_FILENAME,
                    cnwheat_simulation.Simulation.HIDDENZONE_T_INDEXES,
                ),
                (
                    self.elements_all_data_list,
                    self.ELEMENTS_OUTPUTS_FILENAME,
                    cnwheat_simulation.Simulation.ELEMENTS_T_INDEXES,
                ),
                (self.soils_all_data_list, self.SOILS_OUTPUTS_FILENAME, cnwheat_simulation.Simulation.SOILS_T_INDEXES),
            ):
                data_filepath = os.path.join(outputs_brut, outputs_filename)
                outputs_df = pandas.concat(outputs_df_list, keys=self.all_simulation_steps, sort=False)
                outputs_df.reset_index(0, inplace=True)
                outputs_df.rename({"level_0": "t"}, axis=1, inplace=True)
                outputs_df = outputs_df.reindex(
                    index_columns + outputs_df.columns.difference(index_columns).tolist(), axis=1, copy=False
                )
                outputs_df.fillna(value=numpy.nan, inplace=True)  # Convert back None to NaN
                save_df_to_csv(outputs_df, data_filepath, PRECISION)
                outputs_file_basename = outputs_filename.split(".")[0]
                outputs_df_dict[outputs_file_basename] = outputs_df.reset_index()

        # construit un premier postprocessing
        if run_postprocessing:
            POSTPROCESSING_DIRPATH = os.path.join(self.out_folder, "postprocessing")

            time_grid = list(outputs_df_dict.values())[0].t
            delta_t = (time_grid.loc[1] - time_grid.loc[0]) * self.HOUR_TO_SECOND_CONVERSION_FACTOR

            axes_postprocessing_file_basename = self.AXES_POSTPROCESSING_FILENAME.split(".")[0]
            hiddenzones_postprocessing_file_basename = self.HIDDENZONES_POSTPROCESSING_FILENAME.split(".")[0]
            organs_postprocessing_file_basename = self.ORGANS_POSTPROCESSING_FILENAME.split(".")[0]
            elements_postprocessing_file_basename = self.ELEMENTS_POSTPROCESSING_FILENAME.split(".")[0]
            soils_postprocessing_file_basename = self.SOILS_POSTPROCESSING_FILENAME.split(".")[0]

            postprocessing_df_dict = {}
            (
                postprocessing_df_dict[axes_postprocessing_file_basename],
                postprocessing_df_dict[hiddenzones_postprocessing_file_basename],
                postprocessing_df_dict[organs_postprocessing_file_basename],
                postprocessing_df_dict[elements_postprocessing_file_basename],
                postprocessing_df_dict[soils_postprocessing_file_basename],
            ) = cnwheat_facade.CNWheatFacade.postprocessing(
                axes_outputs_df=outputs_df_dict[self.AXES_OUTPUTS_FILENAME.split(".")[0]],
                hiddenzone_outputs_df=outputs_df_dict[self.HIDDENZONES_OUTPUTS_FILENAME.split(".")[0]],
                organs_outputs_df=outputs_df_dict[self.ORGANS_OUTPUTS_FILENAME.split(".")[0]],
                elements_outputs_df=outputs_df_dict[self.ELEMENTS_OUTPUTS_FILENAME.split(".")[0]],
                soils_outputs_df=outputs_df_dict[self.SOILS_OUTPUTS_FILENAME.split(".")[0]],
                delta_t=delta_t,
            )
           

            for postprocessing_file_basename, postprocessing_filename in (
                (axes_postprocessing_file_basename, self.AXES_POSTPROCESSING_FILENAME),
                (hiddenzones_postprocessing_file_basename, self.HIDDENZONES_POSTPROCESSING_FILENAME),
                (organs_postprocessing_file_basename, self.ORGANS_POSTPROCESSING_FILENAME),
                (elements_postprocessing_file_basename, self.ELEMENTS_POSTPROCESSING_FILENAME),
                (soils_postprocessing_file_basename, self.SOILS_POSTPROCESSING_FILENAME),
            ):
                postprocessing_filepath = os.path.join(POSTPROCESSING_DIRPATH, postprocessing_filename)
                postprocessing_df_dict[postprocessing_file_basename].to_csv(
                    postprocessing_filepath, na_rep="NA", index=False, float_format="%.{}f".format(PRECISION)
                )

                if run_graphs:

                    GRAPHS_DIRPATH = os.path.join(self.out_folder, "graphs")
                     # --- Generate graphs from postprocessing files
                    plt.ioff()
                    df_elt = postprocessing_df_dict[elements_postprocessing_file_basename]
                    df_SAM = pandas.read_csv(os.path.join(outputs_brut, self.AXES_OUTPUTS_FILENAME))

                    cnwheat_facade.CNWheatFacade.graphs(axes_postprocessing_df=postprocessing_df_dict[axes_postprocessing_file_basename],
                                                        hiddenzones_postprocessing_df=postprocessing_df_dict[hiddenzones_postprocessing_file_basename],
                                                        organs_postprocessing_df=postprocessing_df_dict[organs_postprocessing_file_basename],
                                                        elements_postprocessing_df=postprocessing_df_dict[elements_postprocessing_file_basename],
                                                        soils_postprocessing_df=postprocessing_df_dict[soils_postprocessing_file_basename],
                                                        meteo_data=self.meteo, graphs_dirpath=GRAPHS_DIRPATH)

                                
                    # --- Additional graphs
                    from cnwheat import tools as cnwheat_tools
                    colors = ['blue', 'darkorange', 'green', 'red', 'darkviolet', 'gold', 'magenta', 'brown', 'darkcyan', 'grey', 'lime']
                    colors = colors + colors

                    # 0) Phyllochron
                    df_SAM = df_SAM[df_SAM['axis'] == 'MS']
                    df_hz = postprocessing_df_dict[hiddenzones_postprocessing_file_basename]
                    grouped_df = df_hz[df_hz['axis'] == 'MS'].groupby(['plant', 'metamer'])[['t', 'leaf_is_emerged']]
                    leaf_emergence = {}
                    for group_name, data in grouped_df:
                        plant, metamer = group_name[0], group_name[1]
                        if metamer == 3 or True not in data['leaf_is_emerged'].unique():
                            continue
                        leaf_emergence_t = data[data['leaf_is_emerged'] == True].iloc[0]['t']
                        leaf_emergence[(plant, metamer)] = leaf_emergence_t

                    phyllochron = {'plant': [], 'metamer': [], 'phyllochron': []}
                    for key, leaf_emergence_t in sorted(leaf_emergence.items()):
                        plant, metamer = key[0], key[1]
                        if metamer == 4:
                            continue
                        phyllochron['plant'].append(plant)
                        phyllochron['metamer'].append(metamer)
                        prev_leaf_emergence_t = leaf_emergence[(plant, metamer - 1)]
                        if df_SAM[(df_SAM['t'] == leaf_emergence_t) | (df_SAM['t'] == prev_leaf_emergence_t)].sum_TT.count() == 2:
                            phyllo_DD = df_SAM[(df_SAM['t'] == leaf_emergence_t)].sum_TT.values[0] - df_SAM[(df_SAM['t'] == prev_leaf_emergence_t)].sum_TT.values[0]
                        else:
                            phyllo_DD = numpy.nan
                        phyllochron['phyllochron'].append(phyllo_DD)

                    if len(phyllochron['metamer']) > 0:
                        fig, ax = plt.subplots()
                        plt.xlim((int(min(phyllochron['metamer']) - 1), int(max(phyllochron['metamer']) + 1)))
                        plt.ylim(ymin=0, ymax=150)
                        ax.plot(phyllochron['metamer'], phyllochron['phyllochron'], color='b', marker='o')
                        for i, j in zip(phyllochron['metamer'], phyllochron['phyllochron']):
                            ax.annotate(str(int(round(j, 0))), xy=(i, j + 2), ha='center')
                        ax.set_xlabel('Leaf number')
                        ax.set_ylabel('Phyllochron (Degree Day)')
                        ax.set_title('phyllochron')
                        plt.savefig(os.path.join(GRAPHS_DIRPATH, 'phyllochron' + '.PNG'))
                        plt.close()

                    # 1) Comparison Dimensions with Ljutovac 2002
                    data_obs = pandas.read_csv(r'inputs_fspmwheat\Ljutovac2002.csv')
                    bchmk = data_obs
                    res = pandas.read_csv(os.path.join(outputs_brut, self.HIDDENZONES_OUTPUTS_FILENAME))
                    res = res[(res['axis'] == 'MS') & (res['plant'] == 1) & ~numpy.isnan(res.leaf_Lmax)].copy()
                    res_IN = res[~ numpy.isnan(res.internode_Lmax)]
                    last_value_idx = res.groupby(['metamer'])['t'].transform(max) == res['t']
                    res = res[last_value_idx].copy()
                    res['lamina_Wmax'] = res.leaf_Wmax
                    res['lamina_W_Lg'] = res.leaf_Wmax / res.lamina_Lmax
                    bchmk = bchmk.loc[bchmk.metamer >= min(res.metamer)]
                    bchmk['lamina_W_Lg'] = bchmk.lamina_Wmax / bchmk.lamina_Lmax
                    last_value_idx = res_IN.groupby(['metamer'])['t'].transform(max) == res_IN['t']
                    res_IN = res_IN[last_value_idx].copy()
                    res = res[['metamer', 'leaf_Lmax', 'lamina_Lmax', 'sheath_Lmax', 'lamina_Wmax', 'lamina_W_Lg', 'SSLW', 'LSSW']].merge(res_IN[['metamer', 'internode_Lmax']], left_on='metamer',
                                                                                                                                        right_on='metamer', how='outer').copy()

                    var_list = ['leaf_Lmax', 'lamina_Lmax', 'sheath_Lmax', 'lamina_Wmax', 'internode_Lmax']
                    for var in list(var_list):
                        fig, ax = plt.subplots()
                        plt.xlim((int(min(res.metamer) - 1), int(max(res.metamer) + 1)))
                        plt.ylim(ymin=0, ymax=numpy.nanmax(list(res[var] * 100 * 1.05) + list(bchmk[var] * 1.05)))

                        tmp = res[['metamer', var]].drop_duplicates()

                        line1 = ax.plot(tmp.metamer, tmp[var] * 100, color='c', marker='o')
                        line2 = ax.plot(bchmk.metamer, bchmk[var], color='orange', marker='o')

                        ax.set_ylabel(var + ' (cm)')
                        ax.set_title(var)
                        ax.legend((line1[0], line2[0]), ('Simulation', 'Ljutovac 2002'), loc=2)
                        plt.savefig(os.path.join(GRAPHS_DIRPATH, var + '.PNG'))
                        plt.close()

                    var = 'lamina_W_Lg'
                    fig, ax = plt.subplots()
                    plt.xlim((int(min(res.metamer) - 1), int(max(res.metamer) + 1)))
                    plt.ylim(ymin=0, ymax=numpy.nanmax(list(res[var] * 1.05) + list(bchmk[var] * 1.05)))
                    tmp = res[['metamer', var]].drop_duplicates()
                    line1 = ax.plot(tmp.metamer, tmp[var], color='c', marker='o')
                    line2 = ax.plot(bchmk.metamer, bchmk[var], color='orange', marker='o')
                    ax.set_ylabel(var)
                    ax.set_title(var)
                    ax.legend((line1[0], line2[0]), ('Simulation', 'Ljutovac 2002'), loc=2)
                    plt.savefig(os.path.join(GRAPHS_DIRPATH, var + '.PNG'))
                    plt.close()

                    # 1bis) Comparison Structural Masses vs. adaptation from Bertheloot 2008

                    # SSLW Laminae
                    bchmk = pandas.DataFrame.from_dict({1: 15, 2: 23, 3: 25, 4: 18, 5: 22, 6: 25, 7: 20, 8: 23, 9: 26, 10: 28, 11: 31}, orient='index').rename(columns={0: 'SSLW'})
                    bchmk.index.name = 'metamer'
                    bchmk = bchmk.reset_index()
                    bchmk = bchmk[bchmk.metamer >= min(res.metamer)]

                    fig, ax = plt.subplots()
                    plt.xlim((int(min(res.metamer) - 1), int(max(res.metamer) + 1)))
                    plt.ylim(ymin=0, ymax=50)

                    tmp = res[['metamer', 'SSLW']].drop_duplicates()

                    line1 = ax.plot(tmp.metamer, tmp.SSLW, color='c', marker='o')
                    line2 = ax.plot(bchmk.metamer, bchmk.SSLW, color='orange', marker='o')

                    ax.set_ylabel('Structural Specific Lamina Weight (g.m-2)')
                    ax.set_title('Structural Specific Lamina Weight')
                    ax.legend((line1[0], line2[0]), ('Simulation', 'adapated from Bertheloot 2008'), loc=3)
                    plt.savefig(os.path.join(GRAPHS_DIRPATH, 'SSLW.PNG'))
                    plt.close()

                    # LWS Sheaths
                    bchmk = pandas.DataFrame.from_dict({1: 0.08, 2: 0.09, 3: 0.11, 4: 0.18, 5: 0.17, 6: 0.21, 7: 0.24, 8: 0.4, 9: 0.5, 10: 0.55, 11: 0.65}, orient='index').rename(columns={0: 'LSSW'})
                    bchmk.index.name = 'metamer'
                    bchmk = bchmk.reset_index()
                    bchmk = bchmk[bchmk.metamer >= min(res.metamer)]

                    fig, ax = plt.subplots()
                    plt.xlim((int(min(res.metamer) - 1), int(max(res.metamer) + 1)))
                    plt.ylim(ymin=0, ymax=0.8)

                    tmp = res[['metamer', 'LSSW']].drop_duplicates()

                    line1 = ax.plot(tmp.metamer, tmp.LSSW, color='c', marker='o')
                    line2 = ax.plot(bchmk.metamer, bchmk.LSSW, color='orange', marker='o')

                    ax.set_ylabel('Lineic Structural Sheath Weight (g.m-1)')
                    ax.set_title('Lineic Structural Sheath Weight')
                    ax.legend((line1[0], line2[0]), ('Simulation', 'adapated from Bertheloot 2008'), loc=2)
                    plt.savefig(os.path.join(GRAPHS_DIRPATH, 'LSSW.PNG'))
                    plt.close()

                    # 2) LAI
                    df_elt['green_area_rep'] = df_elt.green_area * df_elt.nb_replications
                    grouped_df = df_elt[(df_elt.axis == 'MS') & (df_elt.element == 'LeafElement1')].groupby(['t', 'plant'])
                    LAI_dict = {'t': [], 'plant': [], 'LAI': []}
                    for name, data in grouped_df:
                        t, plant = name[0], name[1]
                        LAI_dict['t'].append(t)
                        LAI_dict['plant'].append(plant)
                        LAI_dict['LAI'].append(data['green_area_rep'].sum() * self.plant_density[plant])

                    cnwheat_tools.plot_cnwheat_ouputs(pandas.DataFrame(LAI_dict), 't', 'LAI', x_label='Time (Hour)', y_label='LAI', plot_filepath=os.path.join(GRAPHS_DIRPATH, 'LAI.PNG'), explicit_label=False)

                    # 3) RER during the exponentiel-like phase

                    # - RER parameters
                    rer_param = dict((k, v) for k, v in elongwheat_parameters.RERmax.items())

                    # - Simulated RER

                    # import simulation outputs
                    data_RER = pandas.read_csv(os.path.join(outputs_brut, self.HIDDENZONES_OUTPUTS_FILENAME))
                    data_RER = data_RER[(data_RER.axis == 'MS') & (data_RER.metamer >= 4)].copy()
                    data_RER.sort_values(['t', 'metamer'], inplace=True)
                    data_teq = pandas.read_csv(os.path.join(outputs_brut, self.AXES_OUTPUTS_FILENAME))
                    data_teq = data_teq[data_teq.axis == 'MS'].copy()

                    # - Time previous leaf emergence
                    tmp = data_RER[data_RER.leaf_is_emerged]
                    leaf_em = tmp.groupby('metamer', as_index=False)['t'].min()
                    leaf_em['t_em'] = leaf_em.t
                    prev_leaf_em = leaf_em
                    prev_leaf_em.metamer = leaf_em.metamer + 1

                    data_RER2 = pandas.merge(data_RER, prev_leaf_em[['metamer', 't_em']], on='metamer')
                    data_RER2 = data_RER2[data_RER2.t <= data_RER2.t_em]

                    # - SumTimeEq
                    data_teq['SumTimeEq'] = numpy.cumsum(data_teq.delta_teq)
                    data_RER3 = pandas.merge(data_RER2, data_teq[['t', 'SumTimeEq']], on='t')

                    # - logL
                    data_RER3['logL'] = numpy.log(data_RER3.leaf_L)

                    # - Estimate RER
                    RER_sim = {}
                    for leaf in data_RER3.metamer.drop_duplicates():
                        Y = data_RER3.logL[data_RER3.metamer == leaf]
                        X = data_RER3.SumTimeEq[data_RER3.metamer == leaf]
                        X = sm.add_constant(X)
                        mod = sm.OLS(Y, X)
                        fit_RER = mod.fit()
                        RER_sim[leaf] = fit_RER.params['SumTimeEq']

                    # - Graph
                    fig, ax1 = plt.subplots()
                    ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))

                    x, y = zip(*sorted(RER_sim.items()))
                    ax1.plot(x, y, label=r'Simulated RER', linestyle='-', color='g')
                    ax1.errorbar(data_obs.metamer, data_obs.RER, yerr=data_obs.RER_confint, marker='o', color='g', linestyle='', label="Observed RER", markersize=2)
                    ax1.plot(list(rer_param.keys()), list(rer_param.values()), marker='*', color='k', linestyle='', label="Model parameters")

                    # Formatting
                    ax1.set_ylabel(u'Relative Elongation Rate at 12�C (s$^{-1}$)')
                    ax1.legend(prop={'size': 12}, bbox_to_anchor=(0.05, .6, 0.9, .5), loc='upper center', ncol=3, mode="expand", borderaxespad=0.)
                    ax1.legend(loc='upper left')
                    ax1.set_xlabel('Phytomer rank')
                    ax1.set_ylim(bottom=0., top=6e-6)
                    ax1.set_xlim(left=4)
                    plt.savefig(os.path.join(GRAPHS_DIRPATH, 'RER_comparison.PNG'), format='PNG', bbox_inches='tight', dpi=200)
                    plt.close()

                    # 4) Total C production vs. Root C allcoation
                    df_org = postprocessing_df_dict[organs_postprocessing_file_basename]
                    df_roots = df_org[df_org['organ'] == 'roots'].copy()
                    df_roots['day'] = df_roots['t'] // 24 + 1
                    df_roots['Unloading_Sucrose_tot'] = df_roots['Unloading_Sucrose'] * df_roots['mstruct']
                    Unloading_Sucrose_tot = df_roots.groupby(['day'])['Unloading_Sucrose_tot'].agg('sum')
                    days = df_roots['day'].unique()

                    df_axe = postprocessing_df_dict[axes_postprocessing_file_basename]
                    df_axe['day'] = df_axe['t'] // 24 + 1
                    Total_Photosynthesis = df_axe.groupby(['day'])['Tillers_Photosynthesis'].agg('sum')

                    df_elt = postprocessing_df_dict[elements_postprocessing_file_basename]
                    df_elt['day'] = df_elt['t'] // 24 + 1
                    df_elt['sum_respi_tillers'] = df_elt['sum_respi'] * df_elt['nb_replications']
                    Shoot_respiration = df_elt.groupby(['day'])['sum_respi_tillers'].agg('sum')
                    Net_Photosynthesis = Total_Photosynthesis - Shoot_respiration

                    share_net_roots_live = Unloading_Sucrose_tot / Net_Photosynthesis * 100

                    fig, ax = plt.subplots()
                    line1 = ax.plot(days, Net_Photosynthesis, label=u'Net_Photosynthesis')
                    line2 = ax.plot(days, Unloading_Sucrose_tot, label=u'C unloading to roots')

                    ax2 = ax.twinx()
                    line3 = ax2.plot(days, share_net_roots_live, label=u'Net C Shoot production sent to roots (%)', color='red')

                    lines = line1 + line2 + line3
                    labs = [line.get_label() for line in lines]
                    ax.legend(lines, labs, loc='center left', prop={'size': 10}, framealpha=0.5, bbox_to_anchor=(1, 0.815), borderaxespad=0.)

                    ax.set_xlabel('Days')
                    ax2.set_ylim([0, 200])
                    ax.set_ylabel(u'C (�mol C.day$^{-1}$ )')
                    ax2.set_ylabel(u'Ratio (%)')
                    ax.set_title('C allocation to roots')
                    plt.savefig(os.path.join(GRAPHS_DIRPATH, 'C_allocation.PNG'), dpi=200, format='PNG', bbox_inches='tight')

                    # 5) C usages relatif to Net Photosynthesis
                    df_org = postprocessing_df_dict[organs_postprocessing_file_basename]
                    df_roots = df_org[df_org['organ'] == 'roots'].copy()
                    df_roots['day'] = df_roots['t'] // 24 + 1
                    df_phloem = df_org[df_org['organ'] == 'phloem'].copy()
                    df_phloem['day'] = df_phloem['t'] // 24 + 1

                    AMINO_ACIDS_C_RATIO = 4.15  #: Mean number of mol of C in 1 mol of the major amino acids of plants (Glu, Gln, Ser, Asp, Ala, Gly)
                    AMINO_ACIDS_N_RATIO = 1.25  #: Mean number of mol of N in 1 mol of the major amino acids of plants (Glu, Gln, Ser, Asp, Ala, Gly)

                    # Photosynthesis
                    df_elt['Photosynthesis_tillers'] = df_elt['Photosynthesis'].fillna(0) * df_elt['nb_replications'].fillna(1.)
                    Tillers_Photosynthesis_Ag = df_elt.groupby(['t'], as_index=False).agg({'Photosynthesis_tillers': 'sum'})
                    C_usages = pandas.DataFrame({'t': Tillers_Photosynthesis_Ag['t']})
                    C_usages['C_produced'] = numpy.cumsum(Tillers_Photosynthesis_Ag.Photosynthesis_tillers)
                    C_usages['C_produced'] = numpy.where(C_usages['C_produced']==0, 1, C_usages['C_produced'])

                    # Respiration
                    C_usages['Respi_roots'] = numpy.cumsum(df_axe.C_respired_roots)
                    C_usages['Respi_shoot'] = numpy.cumsum(df_axe.C_respired_shoot)

                    # Exudation
                    C_usages['exudation'] = numpy.cumsum(df_axe.C_exudated.fillna(0))

                    # Structural growth
                    C_consumption_mstruct_roots = df_roots.sucrose_consumption_mstruct.fillna(0) + df_roots.AA_consumption_mstruct.fillna(0) * AMINO_ACIDS_C_RATIO / AMINO_ACIDS_N_RATIO
                    C_usages['Structure_roots'] = numpy.cumsum(C_consumption_mstruct_roots.reset_index(drop=True))

                    df_hz['C_consumption_mstruct'] = df_hz.sucrose_consumption_mstruct.fillna(0) + df_hz.AA_consumption_mstruct.fillna(0) * AMINO_ACIDS_C_RATIO / AMINO_ACIDS_N_RATIO
                    df_hz['C_consumption_mstruct_tillers'] = df_hz['C_consumption_mstruct'] * df_hz['nb_replications']
                    C_consumption_mstruct_shoot = df_hz.groupby(['t'])['C_consumption_mstruct_tillers'].sum()
                    C_usages['Structure_shoot'] = numpy.cumsum(C_consumption_mstruct_shoot.reset_index(drop=True))

                    # Non structural C
                    df_phloem['C_NS'] = df_phloem.sucrose.fillna(0) + df_phloem.amino_acids.fillna(0) * AMINO_ACIDS_C_RATIO / AMINO_ACIDS_N_RATIO
                    C_NS_phloem_init = df_phloem.C_NS - df_phloem.C_NS[0]
                    C_usages['NS_phloem'] = C_NS_phloem_init.reset_index(drop=True)

                    df_elt['C_NS'] = df_elt.sucrose.fillna(0) + df_elt.fructan.fillna(0) + df_elt.starch.fillna(0) + (
                            df_elt.amino_acids.fillna(0) + df_elt.proteins.fillna(0)) * AMINO_ACIDS_C_RATIO / AMINO_ACIDS_N_RATIO
                    df_elt['C_NS_tillers'] = df_elt['C_NS'] * df_elt['nb_replications'].fillna(1.)
                    C_elt = df_elt.groupby(['t']).agg({'C_NS_tillers': 'sum'})

                    df_hz['C_NS'] = df_hz.sucrose.fillna(0) + df_hz.fructan.fillna(0) + (df_hz.amino_acids.fillna(0) + df_hz.proteins.fillna(0)) * AMINO_ACIDS_C_RATIO / AMINO_ACIDS_N_RATIO
                    df_hz['C_NS_tillers'] = df_hz['C_NS'] * df_hz['nb_replications'].fillna(1.)
                    C_hz = df_hz.groupby(['t']).agg({'C_NS_tillers': 'sum'})

                    df_roots['C_NS'] = df_roots.sucrose.fillna(0) + df_roots.amino_acids.fillna(0) * AMINO_ACIDS_C_RATIO / AMINO_ACIDS_N_RATIO

                    C_NS_autre = df_roots.C_NS.reset_index(drop=True) + C_elt.C_NS_tillers + C_hz.C_NS_tillers
                    C_NS_autre_init = C_NS_autre - C_NS_autre[0]
                    C_usages['NS_other'] = C_NS_autre_init.reset_index(drop=True)

                    # Total
                    C_usages['C_budget'] = (C_usages.Respi_roots + C_usages.Respi_shoot + C_usages.exudation + C_usages.Structure_roots + C_usages.Structure_shoot + C_usages.NS_phloem + C_usages.NS_other) / \
                                        C_usages.C_produced

                    # ----- Graph
                    fig, ax = plt.subplots()
                    ax.plot(C_usages.t, C_usages.Structure_shoot / C_usages.C_produced * 100,
                            label=u'Structural mass - Shoot', color='g')
                    ax.plot(C_usages.t, C_usages.Structure_roots / C_usages.C_produced * 100,
                            label=u'Structural mass - Roots', color='r')
                    ax.plot(C_usages.t, (C_usages.NS_phloem + C_usages.NS_other) / C_usages.C_produced * 100, label=u'Non-structural C', color='darkorange')
                    ax.plot(C_usages.t, (C_usages.Respi_roots + C_usages.Respi_shoot) / C_usages.C_produced * 100, label=u'C loss by respiration', color='b')
                    ax.plot(C_usages.t, C_usages.exudation / C_usages.C_produced * 100, label=u'C loss by exudation', color='c')

                    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
                    ax.set_xlabel('Time (h)')
                    ax.set_ylabel(u'Carbon usages : Photosynthesis (%)')
                    ax.set_ylim(bottom=0, top=100.)

                    fig.suptitle(u'Total cumulated usages are ' + str(round(C_usages.C_budget.tail(1) * 100, 0)) + u' % of Photosynthesis')

                    plt.savefig(os.path.join(GRAPHS_DIRPATH, 'C_usages_cumulated.PNG'), format='PNG', bbox_inches='tight')
                    plt.close()

                    # 6) RUE
                    df_elt['PARa_MJ'] = df_elt['PARa'] * df_elt['green_area'] * df_elt['nb_replications'] * 3600 / 4.6 * 10 ** -6  # Il faudrait idealement utiliser les calculcs green_area et PARa des talles
                    df_elt['RGa_MJ'] = df_elt['PARa'] * df_elt['green_area'] * df_elt['nb_replications'] * 3600 / 2.02 * 10 ** -6  # Il faudrait idealement utiliser les calculcs green_area et PARa des talles
                    PARa = df_elt.groupby(['day'])['PARa_MJ'].agg('sum')
                    PARa_cum = numpy.cumsum(PARa)
                    days = df_elt['day'].unique()

                    sum_dry_mass_shoot = df_axe.groupby(['day'])['sum_dry_mass_shoot'].agg('max').astype('float64')
                    sum_dry_mass = df_axe.groupby(['day'])['sum_dry_mass'].agg('max').astype('float64')

                    RUE_shoot = numpy.polyfit(PARa_cum, sum_dry_mass_shoot, 1)[0]
                    RUE_plant = numpy.polyfit(PARa_cum, sum_dry_mass, 1)[0]

                    fig, ax = plt.subplots()
                    ax.plot(PARa_cum, sum_dry_mass_shoot, label='Shoot dry mass (g)')
                    ax.plot(PARa_cum, sum_dry_mass, label='Plant dry mass (g)')
                    ax.legend(prop={'size': 10}, framealpha=0.5, loc='center left', bbox_to_anchor=(1, 0.815), borderaxespad=0.)
                    ax.set_xlabel('Cumulative absorbed PAR (MJ)')
                    ax.set_ylabel('Dry mass (g)')
                    ax.set_title('RUE')
                    plt.text(max(PARa_cum) * 0.02, max(sum_dry_mass) * 0.95, 'RUE shoot : {0:.2f} , RUE plant : {1:.2f}'.format(round(RUE_shoot, 2), round(RUE_plant, 2)))
                    plt.savefig(os.path.join(GRAPHS_DIRPATH, 'RUE.PNG'), dpi=200, format='PNG', bbox_inches='tight')

                    fig, ax = plt.subplots()
                    ax.plot(days, sum_dry_mass_shoot, label='Shoot dry mass (g)')
                    ax.plot(days, sum_dry_mass, label='Plant dry mass (g)')
                    ax.plot(days, PARa_cum, label='Cumulative absorbed PAR (MJ)')
                    ax.legend(prop={'size': 10}, framealpha=0.5, loc='center left', bbox_to_anchor=(1, 0.815), borderaxespad=0.)
                    ax.set_xlabel('Days')
                    ax.set_title('RUE investigations')
                    plt.savefig(os.path.join(GRAPHS_DIRPATH, 'RUE2.PNG'), dpi=200, format='PNG', bbox_inches='tight')

                    # 7) Sum thermal time
                    df_SAM = df_SAM[df_SAM['axis'] == 'MS']
                    fig, ax = plt.subplots()
                    ax.plot(df_SAM['t'], df_SAM['sum_TT'])
                    ax.set_xlabel('Hours')
                    ax.set_ylabel('Thermal Time')
                    ax.set_title('Thermal Time')
                    plt.savefig(os.path.join(GRAPHS_DIRPATH, 'SumTT.PNG'), dpi=200, format='PNG', bbox_inches='tight')

                    # 7) Residual N : ratio_N_mstruct_max
                    df_elt_outputs = pandas.read_csv(os.path.join(outputs_brut, self.ELEMENTS_OUTPUTS_FILENAME))
                    df_elt_outputs = df_elt_outputs.loc[df_elt_outputs.axis == 'MS']
                    df_elt_outputs = df_elt_outputs.loc[df_elt_outputs.mstruct != 0]
                    df_elt_outputs['N_content_total'] = df_elt_outputs['N_content_total'] * 100
                    x_name = 't'
                    x_label = 'Time (Hour)'
                    graph_variables_ph_elements = {'N_content_total': u'N content in green + senesced tissues (% mstruct)'}
                    for org_ph in (['blade'], ['sheath'], ['internode'], ['peduncle', 'ear']):
                        for variable_name, variable_label in graph_variables_ph_elements.items():
                            graph_name = variable_name + '_' + '_'.join(org_ph) + '.PNG'
                            cnwheat_tools.plot_cnwheat_ouputs(df_elt_outputs,
                                                            x_name=x_name,
                                                            y_name=variable_name,
                                                            x_label=x_label,
                                                            y_label=variable_label,
                                                            colors=[colors[i - 1] for i in df_elt_outputs.metamer.unique().tolist()],
                                                            filters={'organ': org_ph},
                                                            plot_filepath=os.path.join(GRAPHS_DIRPATH, graph_name),
                                                            explicit_label=False)

        self.outputs_df_dict = outputs_df_dict

    def update_Nitrates_cnwheat_mtg(self):
        """Transfers N uptake result in model MTG
        """        
        mtg_plants_iterator = self.g.components_iter(self.g.root)
        for plant in self.cnwheat_facade_.population.plants:
            cnwheat_plant_index = plant.index
            while True:
                mtg_plant_vid = next(mtg_plants_iterator)
                if int(self.g.index(mtg_plant_vid)) == cnwheat_plant_index:
                    break
            mtg_axes_iterator = self.g.components_iter(mtg_plant_vid)
            for axis in plant.axes:
                # update in cnwheat_facade_
                axis.roots.__dict__["Uptake_Nitrates"] = self.uptake_nitrate_hour

                # update Nitrates uptake in MTG
                cnwheat_axis_label = axis.label
                while True:
                    mtg_axis_vid = next(mtg_axes_iterator)
                    if self.g.label(mtg_axis_vid) == cnwheat_axis_label:
                        break
                mtg_roots_properties = self.g.get_vertex_property(mtg_axis_vid)["roots"]
                mtg_roots_properties["Uptake_Nitrates"] = self.uptake_nitrate_hour

    def compute_N_content_roots(self):
        """Get and compute N content roots at current timestep

        Returns
        -------
        float
            N content in roots [0-1]
        """        
        organs_df = self.organs_all_data_list[-1]
        hiddenzones_df = self.hiddenzones_all_data_list[-1]
        elements_df = self.elements_all_data_list[-1]

        # constantes
        C_MOLAR_MASS = cnwheat_model.EcophysiologicalConstants.C_MOLAR_MASS
        N_MOLAR_MASS = cnwheat_model.EcophysiologicalConstants.N_MOLAR_MASS
        HEXOSE_MOLAR_MASS_C_RATIO = cnwheat_model.EcophysiologicalConstants.HEXOSE_MOLAR_MASS_C_RATIO
        NITRATES_MOLAR_MASS_N_RATIO = cnwheat_model.EcophysiologicalConstants.NITRATES_MOLAR_MASS_N_RATIO
        AMINO_ACIDS_MOLAR_MASS_N_RATIO = cnwheat_model.EcophysiologicalConstants.AMINO_ACIDS_MOLAR_MASS_N_RATIO

        # récupère grandeurs
        structure = organs_df.fillna(0)["structure"]
        sucrose = organs_df.fillna(0)["sucrose"]
        starch = organs_df.fillna(0)["starch"]
        nitrates = organs_df.fillna(0)["nitrates"]
        amino_acids = organs_df.fillna(0)["amino_acids"]
        proteins = organs_df.fillna(0)["proteins"]
        mstruct = organs_df.fillna(0)["mstruct"]
        Nstruct = organs_df.fillna(0)["Nstruct"]

        ## DRY MASS
        organs_df["sum_dry_mass"] = (
            ((structure + starch) * 1e-6 * C_MOLAR_MASS / HEXOSE_MOLAR_MASS_C_RATIO)
            + (sucrose * 1e-6 * C_MOLAR_MASS) / HEXOSE_MOLAR_MASS_C_RATIO
            + (starch * 1e-6 * C_MOLAR_MASS) / HEXOSE_MOLAR_MASS_C_RATIO
            + (nitrates * 1e-6 * N_MOLAR_MASS) / NITRATES_MOLAR_MASS_N_RATIO
            + (amino_acids * 1e-6 * N_MOLAR_MASS) / AMINO_ACIDS_MOLAR_MASS_N_RATIO
            + (proteins * 1e-6 * N_MOLAR_MASS) / AMINO_ACIDS_MOLAR_MASS_N_RATIO
            + mstruct
        )
        dry_mass_roots = (
            organs_df[(organs_df["organ"] == "roots")].groupby(["plant", "axis"])["sum_dry_mass"].agg("sum")
        )

        # calcule de la proportion partie roots du phloem
        hz_df_MS = hiddenzones_df[hiddenzones_df["axis"] == "MS"].copy()
        elt_df_MS = elements_df[elements_df["axis"] == "MS"].copy()
        hz_df_MS["mstruct_tillers"] = hz_df_MS["mstruct"] * hz_df_MS["nb_replications"]
        elt_df_MS["mstruct_tillers"] = elt_df_MS["mstruct"] * elt_df_MS["nb_replications"]
        sum_mstruct_shoot = hz_df_MS.groupby(["plant", "axis"])["mstruct_tillers"].agg("sum") + elt_df_MS.groupby(
            ["plant", "axis"]
        )["mstruct_tillers"].agg("sum")
        sum_mstruct_shoot.fillna(0, inplace=True)
        sum_mstruct_roots = organs_df[(organs_df["organ"] == "roots")].groupby(["plant", "axis"])["mstruct"].agg("sum")
        shoot_roots_mstruct_ratio = sum_mstruct_shoot / sum_mstruct_roots
        phloem_shoot_root = 1 / (1 + 1 / shoot_roots_mstruct_ratio)

        # dry mass de la partie roots du phloem
        sum_dry_mass_phloem = (
            organs_df[(organs_df["organ"] == "phloem")].groupby(["plant", "axis"])["sum_dry_mass"].agg("sum")
        )
        sum_dry_mass_phloem_roots = sum_dry_mass_phloem * (1 - phloem_shoot_root)

        # dry total des roots
        sum_dry_mass_roots = sum_dry_mass_phloem_roots + dry_mass_roots

        ## N MASS
        organs_df["N_g"] = (
            (nitrates * 1e-6 * N_MOLAR_MASS)
            + (amino_acids * 1e-6 * N_MOLAR_MASS)
            + (proteins * 1e-6 * N_MOLAR_MASS)
            + Nstruct
        )

        # N mass du phloem roots
        sum_N_g_phloem = organs_df[(organs_df["organ"] == "phloem")].groupby(["plant", "axis"])["N_g"].agg("sum")
        sum_N_g_phloem_shoot = sum_N_g_phloem * (1 - phloem_shoot_root)

        # N mass du reste des roots
        sum_N_g_roots = organs_df[(organs_df["organ"] == "roots")].groupby(["plant", "axis"])["N_g"].agg("sum")

        # N mass totale des roots
        sum_N_g_roots = sum_N_g_roots + sum_N_g_phloem_shoot

        # finalement pourcentage de N dans les racines
        N_content_roots = sum_N_g_roots / sum_dry_mass_roots * 100

        return N_content_roots.values[0]

    def compute_roots_length(self, soil_wrapper, planter:Planter, t=None):
        """Compute roots length from roots mass

        Note
        ----
        soil_dimensions : [z, x, y]

        Parameters
        ----------
        soil_wrapper : Soil_wrapper
            for soil3ds dimensions and get plant positions in soil voxels grid
        planter : Planter
            for plant positions

        Returns
        -------
        numpy.array
            4-dimensions array of roots length per plant per soil voxel
        """        

        # une seule plante dans cnwheat
        roots_mass = [0]
        for i, plant in enumerate(self.cnwheat_facade_.population.plants):
            for axis in plant.axes:
                roots_mass[i] += axis.roots.mstruct  # masse en g

        positions = planter.wheat_positions[self.wheat_index]
        self.compute_SRL_wheat(roots_mass[0])

        # longueur spécifique x masse en gramme/nbplantes
        ls_roots = []
        for p in positions:
            # on répartit de manière homogène les racines à travers les couches du sol
            # convertit m en cm # --> peut etre en metre finalement
            ix, iy = soil_wrapper.whichvoxel_xy(p)
            roots_length_per_voxel = self.rootsdistribution(roots_mass[0], ix, iy, soil_wrapper, self.rootdistribtype) 
            ls_roots.append(roots_length_per_voxel)

        return ls_roots

    def rootsdistribution(self, roots_mass, ix, iy, soil_wrapper, rootdistribtype="homogeneous", t=None):
        """Distributes roots length in soil voxels.

        Note
        ----
        Currently only homogeneous repartition from ground to deep layers

        Parameters
        ----------
        roots_mass : float
            roots mass in g
        ix : int
            voxel index on x axis
        iy : int
            voxel index on y axis
        soil_wrapper : Soil_wrapper
            for soil voxels dimensions
        rootdistribtype : str, optional
            for futur improvments and other distribution types, by default "homogeneous"

        Returns
        -------
        numpy.array
            3-dimensions array of roots length in each soil voxel
        """        
        roots_length_per_voxel = numpy.zeros(soil_wrapper.soil_dimensions)
        if rootdistribtype == "homogeneous":
            roots_length_per_voxel[:, :, :] = (roots_mass * self.SRL) / (soil_wrapper.soil_dimensions[0]*soil_wrapper.soil_dimensions[1]*soil_wrapper.soil_dimensions[2])

        if rootdistribtype == "column":
            roots_length_per_voxel[:, ix, iy] = (roots_mass * self.SRL) / soil_wrapper.soil_dimensions[0]
        
        if rootdistribtype == "bound":
            for iz in range (self.roots_bound-1):
                roots_length_per_voxel[iz, :, :] = (roots_mass * self.SRL) / (self.roots_bound*soil_wrapper.soil_dimensions[1]*soil_wrapper.soil_dimensions[2])

        if rootdistribtype == "profile":
            #pas fonctionnel, en l'état dilue les racines jusqu'à profondeur max fixe de 150cm
            self.compute_root_profile(roots_mass,soil_wrapper) #computes self.root_profile, 1D array of length soil_wrapper.soil_dimensions[0] (z) containing the root mass per layer
         
            for iz in range(soil_wrapper.soil_dimensions[0]):
                roots_length_per_voxel[iz, :, :] = self.root_profile[iz] * self.SRL / (soil_wrapper.soil_dimensions[1]*soil_wrapper.soil_dimensions[2])
        
        return roots_length_per_voxel

    def compute_SRL_wheat(self, mass_roots):
        """Dynamic specific root length according to roots mass

        Here, SRL is linear following roots mass

        Parameters
        ----------
        mass_roots : float
            roots mass in g
        """        
        if mass_roots < 0.6:
            a = 334
            self.SRL =  a * mass_roots + 60
            # a = 8.247933150630281 # math.log(141)/0.6
            # return numpy.exp(mass_roots * a) + 59
        else:
            self.SRL = 200

    def compute_root_profile(self, roots_mass, soil_wrapper):
        """Dynamic rooting depth according to roots mass

        root profile from mass using eq from Fan et al. 2016 

        Parameters
        ----------
        mass_roots : float
            roots mass in g
        """
        da = 17.2
        c = -1.286 
        dmax = 150.4

        #dmax = 1 * day +20 

        self.root_profile = numpy.zeros(soil_wrapper.soil_dimensions[0])

        for iz in range(len(self.root_profile)):
            d= (iz+1)*soil_wrapper.soil.dxyz[2][0] # depth in cm (iz+1 because iz starts at 0)
            self.root_profile[iz] = (1 / (1 + (d/da)**c) + (1 - 1/(1+ (dmax/da)**c))*(d/dmax) ) * roots_mass 
            # Calculate the root mass at depth d for a total root mass of roots_mass

            # At this stage, root_profile is cumulative; we need to convert it to incremental for root distribution in soil layers

            # Convert the cumulative list to incremental => root mass in each soil layer
        for iz in range(len(self.root_profile) - 1, 0, -1):
            self.root_profile[iz] = self.root_profile[iz] - self.root_profile[iz - 1]
         

    def compute_plants_light_interception(self, plant_leaf_area, soil_energy):
        """Computes light capacity interception for each plant

        Parameters
        ----------
        plant_leaf_area : float
            leaf area of one plant (each plant has the same leaf area as they come from the same mean geometric plant)
        soil_energy : float
            soil light interception [0-1]

        Returns
        -------
        list of float
            light capacity interception for each plant
        """        
        # conversion
        c = (3600 * 24) / 1000000
        # portion au sol pour chaque plante
        return [(1 - soil_energy) * ((plant_leaf_area * c) / self.nb_plants)] * self.nb_plants

    @staticmethod
    def convert_uptake_nitrate(uptake_nitrate):
        """Unit conversion of N uptake from kg to µmol N

        Parameters
        ----------
        uptake_nitrate : float
            N uptake in kg per day

        Returns
        -------
        float
            N uptake in µmol N per hour
        """        
        # masse atomic de l'azote 14.0067 g.mol-1
        atomic_mass_N = 14.0067

        # conversion kg en µmol d'azote
        uptake_nitrate_g = uptake_nitrate * 10**3
        uptake_mol = uptake_nitrate_g / atomic_mass_N
        uptake_micromol = uptake_mol * 10**6

        # daily_uptake_nitrate = 10**9 * uptake_nitrate / atomic_mass_N

        # on répartie la grandeur sur chaque heure
        daily_uptake_nitrate = uptake_micromol / 24

        return daily_uptake_nitrate

    def force_nitrates_uptake(self, t):
        """Transfers forced input N uptake to wheat MTG at timestep t

        Parameters
        ----------
        t : int
            timestep
        """        
        mtg_plants_iterator = self.g.components_iter(self.g.root)
        for plant in self.cnwheat_facade_.population.plants:
            cnwheat_plant_index = plant.index
            while True:
                mtg_plant_vid = next(mtg_plants_iterator)
                if int(self.g.index(mtg_plant_vid)) == cnwheat_plant_index:
                    break
            mtg_axes_iterator = self.g.components_iter(mtg_plant_vid)
            for axis in plant.axes:
                # update Nitrates uptake in dataframe
                group = self.nitrates_uptake_data_grouped.get_group((t, plant.index, axis.label, "roots"))
                nitrates_uptake_data_to_use = (
                    group.loc[group.first_valid_index(), group.columns.intersection(["Uptake_Nitrates"])]
                    .dropna()
                    .to_dict()
                )
                axis.roots.__dict__.update(nitrates_uptake_data_to_use)

                # update Nitrates uptake in MTG
                cnwheat_axis_label = axis.label
                while True:
                    mtg_axis_vid = next(mtg_axes_iterator)
                    if self.g.label(mtg_axis_vid) == cnwheat_axis_label:
                        break
                mtg_roots_properties = self.g.get_vertex_property(mtg_axis_vid)["roots"]
                mtg_roots_properties.update(nitrates_uptake_data_to_use)

    def energy(self, t):
        """Returns meteo radiation at timestep t

        Parameters
        ----------
        t : int
            timestep

        Returns
        -------
        float
            meteo radiation in micromol.m-2.s-1
        """      

        return self.meteo.loc[t, ["PARi"]].iloc[0]

    def doy(self, t, soil3ds=False):
        """Day of the year

        Parameters
        ----------
        t : int
            timestep
        soil3ds : bool, optional
            doy in soil3ds meteo file doesn't reset each year, it keeps increasing, by default False

        Returns
        -------
        int
            day of the year
        """        

        if soil3ds:
            if t > 0:
                if self.meteo.loc[max(0,t - 24), ["DOY"]].iloc[0] > self.meteo.loc[t, ["DOY"]].iloc[0]:
                    self.last_year_doy += self.meteo.loc[t - 24, ["DOY"]].iloc[0]
            return self.meteo.loc[t, ["DOY"]].iloc[0] + self.last_year_doy

        else:
            return self.meteo.loc[t, ["DOY"]].iloc[0]

    def hour(self, t):
        """Hour at timestep t

        Parameters
        ----------
        t : int
            timestep

        Returns
        -------
        int
            hour
        """ 

        return self.meteo.loc[t, ["hour"]].iloc[0]

    def PARi_next_hours(self, t):
        """Meteo PARi at t+1

        Parameters
        ----------
        t : int
            timestep

        Returns
        -------
        float
            meteo PAR at t+1
        """      

        return self.meteo.loc[range(t, t + self.LIGHT_TIMESTEP), ["PARi"]].sum().values[0]

    def next_day_next_hour(self, t):
        """Day of the year at t+1. Used to check if we are at last iteration of the day

        Parameters
        ----------
        t : int
            timestep

        Returns
        -------
        int
            day of the year
        """     

        return self.meteo.loc[t + self.SENESCWHEAT_TIMESTEP, ["DOY"]].iloc[0]


def passive_lighting(data, t, DOY, scene, lighting_wrapper, stems=None):
    """Run the lighting computation step without saving

    Parameters
    ----------
    data : dict
        WheatFspm lighting step results saved
    t : int
        timestep
    DOY : int
        day of the year
    scene : plantgl.Scene
        canopy plantgl scene
    lighting_wrapper : Light_wrapper
        which light model to run
    stems : list, optional
        lits of stems, (specy id, organ id), by default None
    """    

    lighting_wrapper.run(scenes=[scene], day=DOY, parunit="micromol.m-2.s-1", stems=stems)

    results = lighting_wrapper.results_organs()

    para = results["Organ"] * results["Area"]
    para *= 1 / results["Area"].sum()

    data["PARa"].append(para)
    data["t"].append(t)


def import_fertilization(file_name, sheet_name, column_name,wheather_file):
    # Charger le fichier mgmt de l-egume
    legume_mgmt = pandas.read_excel(file_name, sheet_name=sheet_name)

    # Charger le fichier meteo de CNwheat
    df2 = pandas.read_csv(wheather_file)

      # Garder uniquement la première occurrence de chaque DOY dans df2
    df2 = df2.drop_duplicates(subset='DOY', keep='first')

    # Fusionner les deux fichiers Excel en utilisant la colonne DOY comme clé
    df = pandas.merge(legume_mgmt, df2, on='DOY')

    # Créer un dictionnaire à partir des colonnes t et column_name (NO3 ou NH4) du fichier fusionné
    dico = df.set_index('t')[column_name].to_dict()

    # Supprimer les valeurs nulles du dictionnaire
    dico = {k: v for k, v in dico.items() if pandas.notna(v) and v != 0}

     # Diviser tous les éléments du dictionnaire par 140*10^-6
    for key in dico:
        dico[key] /= 140*10**-6

    return dico

