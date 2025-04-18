"""

    contains Planter class

"""

import numpy
import math
import pandas

import openalea.plantgl.all as plantgl
import random

from plantfusion.indexer import Indexer


class Planter:
    """Tool for handling the following in the simulation:
        
        * plant positions
        
        * number of plants 
        
        * soil domain

    You can choose between 3 generation types:

        * "default": each instance of fspm wrapper manages its plant positions, manages only soil domain

        * "random": generates each plant in random positions in a squared soil

        * "row": generates two rows for each plant specy in a squared soil

    Parameters
    ----------
    generation_type : str, optional
        choose between "default", "random" or "row", by default "default"
    indexer : Indexer, optional
        indexer for listing FSPM in the simulation, by default Indexer()
    legume_cote : dict, optional
        precise the length of a soil side in l-egume instances. An entry is {wrapper.name : length in cm}, by default {}
    legume_number_of_plants : dict, optional
        number of plants for each l-egume instances. An entry is {wrapper.name : number of plants}, by default {}
    inter_rows : float, optional
        length between two rows in m, by default 0.15
    plant_density : dict, optional
        number plants in 1 m^2. An entry is {wrapper.name : nb of plants/m² }, by default {1: 250}
    xy_plane : tuple of tuple, optional
        Forced input dimensions of soil domain, ((xmin, ymin), (xmax, ymax)), by default None
    xy_square_length : float, optional
        Side length of the soil domain for "random" generation type, by default 0.5
    translate : dict, optional
        Possibility to translate some of the fspm geometric scenes by a 3d vector. An entry is { wrapper.name : (tx, ty, tz) }, by default None
    noise_plant_positions : float, optional
        noise around generated positions in m, by default 0.0
    save_wheat_positions : bool, optional
        avoid to regenerate wheat positions at each timestep, by default False
    seed : int, optional
        seed for random and numpy, by default None

    """    
    def __init__(
        self,
        generation_type="default",
        indexer=Indexer(),
        legume_cote={},
        legume_number_of_plants={},
        inter_rows=0.15,
        plant_density={1: 250},
        xy_plane=None,
        xy_square_length=0.5,
        translate=None,
        noise_plant_positions=0.0,
        save_wheat_positions=False,
        seed=None,
        positions_grid={},
        n_rows=None, 
        n_cols=None, 
        cell_size=None, 
        col_pattern=None
    ) -> None:
        """Constructor, computes a global soil domain for the simulation

        """        
        self.generation_type = generation_type
        self.plant_density = plant_density
        self.save_wheat_positions = save_wheat_positions
        self.noise_plant_positions = noise_plant_positions
        self.positions_grid = positions_grid
        self.indexer = indexer
        self.number_of_plants: list = [0 for i in indexer.global_order]
        self.domain=None

        # les lsystem l-egume sont par défaut en cm et le reste en m
        self.transformations: dict = {"scenes unit": {}}


        for i in range(len(self.indexer.global_order)):
            self.transformations["scenes unit"][i] = "m"
        for i in self.indexer.legume_index:
            self.transformations["scenes unit"][i] = "cm"

        if generation_type == "default":
            self.__default_preconfigured(legume_cote, inter_rows, plant_density, xy_plane, translate, seed)
            for name, nb_plt in legume_number_of_plants.items():
                self.number_of_plants[indexer.global_order.index(name)] = nb_plt

        elif generation_type == "random":
            self.__random(plant_density, xy_square_length)
            self.type_domain = "mix"

        elif generation_type == "row":
            self.__row(plant_density, inter_rows)
            self.type_domain = "mix"
        
        elif generation_type == "forced":
            self.__forced(positions_grid)
            self.type_domain = "mix"

        elif generation_type == 'row_forced':
            self.domain, positions_grid=self.row_pattern(plant_density, n_rows, inter_rows)
            self.__forced(positions_grid)
            self.type_domain = 'mix'

        elif generation_type == 'grid_forced':
            self.domain, positions_grid=self.grid_pattern(n_cols, n_rows, cell_size, col_pattern , noise=None) #n_cols, n_rows, cell_size, pattern , noise=None
            self.__forced(positions_grid)
            self.type_domain = 'mix'


    def __random(self, plant_density, xy_square_length):
        """Parameters for random generation type

        Note
        ----
        - l-egume: must have 64 plants so it rescale the soil size according to its plant density. Rewrite the cote, nbcote, typearrangement and optdamier options
        - other fspm: compute number of plants according to its plant density

        Parameters
        ----------
        plant_density : dict
            each entry is { wrapper.name : number of plants/m²}
        xy_square_length : float
            soil side length in m
        """        

        self.legume_nbcote = []
        self.wheat_positions = [[] for i in range(len(self.indexer.wheat_names))]
        self.other_positions = [[] for i in range(len(self.indexer.other_names))]

        self.domain = ((0.0, 0.0), (xy_square_length, xy_square_length))
        
        for name, density in plant_density.items():
            if name in self.indexer.legume_names:
                # on réajuste le domaine pour avoir 64 plantes
                xy_square_length = math.sqrt(64 / density)
                self.legume_nbcote.append(8)
                self.domain = ((0.0, 0.0), (xy_square_length, xy_square_length))

                self.legume_typearrangement = "random8"
                # conversion m en cm
                self.legume_cote = xy_square_length * 100

                self.legume_optdamier = 8

                self.number_of_plants[self.indexer.global_order.index(name)] = 64

            else:
                self.number_of_plants[self.indexer.global_order.index(name)] = int(
                    xy_square_length * xy_square_length * density
                )

    def __row(self, plant_density, inter_rows):
        """Parameters for row generation type

        Note
        ----
        - l-egume: use "row4_sp1" arrangement. Rewrite the cote, nbcote, typearrangement and optdamier options
        - other fspm: compute number of plants according to its plant density

        Parameters
        ----------
        plant_density : dict
            each entry is { wrapper.name : number of plants/m²}
        inter_rows : float
            length between two rows in m
        """        
        self.total_n_rows = 2 * len(self.indexer.global_order)

        xy_square_length = inter_rows * self.total_n_rows
        self.domain = ((0.0, 0.0), (xy_square_length, xy_square_length))

        self.inter_rows = inter_rows
        self.legume_nbcote = []
        self.wheat_positions = [[] for i in range(len(self.indexer.wheat_names))]
        self.other_positions = [[] for i in range(len(self.indexer.other_names))]

        for name, density in plant_density.items():

            #nb of plants per fspm instance, based on the sowing density and the size of the domain
            self.number_of_plants[self.indexer.global_order.index(name)] = int(
                xy_square_length * xy_square_length * density
                )

        # each rows are generated on the same positions, we translate each plant specy scenes to create the final scene
        self.transformations["translate"] = {}
        if self.total_n_rows > 4:
            for i in self.indexer.legume_index:
                self.transformations["translate"][i] = (0.0, (i) * inter_rows, 0.0)
            for i in self.indexer.wheat_index:
                self.transformations["translate"][i] = (0.0, (i - 0.5) * inter_rows, 0.0)
            for i in self.indexer.other_index:
                self.transformations["translate"][i] = (0.0, (i - 0.5) * inter_rows, 0.0)

        # only two species
        else:
            # 2 wheats ou combinaison avec 1 ou 2 autres FSPM
            if len(self.indexer.wheat_names) >= 1  or len(self.indexer.other_names) > 1:
                self.transformations["translate"][0] = (0.0, -inter_rows, 0.0)

            # 2 legume
            elif len(self.indexer.legume_names) > 1:
                self.transformations["translate"][1] = (0.0, inter_rows, 0.0)


    def __forced(self,positions_grid):
        """Parameters to force a preset carto into a legume type instance

    
        Parameters
        ----------
        positions : list of arrays
        """        

        #set the size of the domain
        if self.domain is None :

            x_coords=[x for sublist in positions_grid.values() for x, y, z in sublist]
            y_coords=[y for sublist in positions_grid.values() for x, y, z in sublist]

            step_x = round((numpy.unique(x_coords)[1]-numpy.unique(x_coords)[0]),1)
            step_y = round((numpy.unique(y_coords)[1]-numpy.unique(y_coords)[0]),1)

            self.domain = ((0.0, 0.0), (max(x_coords)+step_x/2, max(y_coords)+step_y/2))



        #on crée [[..],[..],..[..]] avec autant de [..] que d'instance de fspm de chaque type 
        self.legume_positions = [[] for i in range(len(self.indexer.legume_names))]
        self.wheat_positions =  [[] for i in range(len(self.indexer.wheat_names))]
        self.other_positions =  [[] for i in range(len(self.indexer.other_names))]

        
        for name, positions in positions_grid.items():
            if name in self.indexer.global_order:
                self.number_of_plants[self.indexer.global_order.index(name)] = len(positions)

                if name in self.indexer.legume_names:
                    i = self.indexer.legume_names.index(name)
                    self.legume_positions[i] = positions
                elif name in self.indexer.wheat_names:
                    i = self.indexer.wheat_names.index(name)
                    self.wheat_positions[i] = positions
                else:
                    i = self.indexer.other_names.index(name)
                    self.other_positions[i] = positions


        #legume parameters are initialised with default values that will get rewritten in legume-wrapper. they are mandatory for L_egume to initalize.
        self.legume_typearrangement = "random8"
        self.legume_nbcote=[8]
        self.legume_cote = (self.domain[1][0]-self.domain[0][0]) *100 #in cm
        self.legume_optdamier = 8

        
    def row_pattern(self, densities, n_rows, inter_rows, offset=None, noise=None):
        # Créer un dictionnaire pour stocker les grilles de coordonnées de chaque catégorie de points
        grids = {}


        # Calculer la largeur du domaine en fonction de l'espacement entre les colonnes le plus élevé et du nombre de colonnes le plus élevé
        # Calculer la largeur totale pour chaque catégorie
        tot_width = {}
        for name in n_rows.keys():
            tot_width[name] = n_rows[name] * inter_rows[name] + (offset[name] if offset is not None and name in offset else inter_rows[name]/2)

        # Trouver la catégorie ayant la plus grande largeur totale
        max_cat = max(tot_width.items(), key=lambda x: x[1])[0]


        max_inter_row = inter_rows[max_cat]
        max_n_rows = n_rows[max_cat]
        common_offset= offset if offset is not None else max_inter_row/(2 * len(self.indexer.global_order) )
        # ^ TEMPORAiRE : suppose que IR est le même pour toutes les cultures et que le nombre de rang est le même pour toutes les cultures
        width = (max_n_rows * max_inter_row)
        # Calculer les coordonnées min et max du domaine carré
        x_min, y_min = 0, 0
        x_max, y_max = width, width

        # Parcourir les catégories de points
        #pour offset automatique si on ne précise pas : le but est de décaller d'1 inter_rang à chaque nouvelle culture 
        rank=0 #va de 0 à nbr de cultures
        
        for name, density in densities.items():

            # Calculer le nombre de plantes dans le domaine en fonction de la densité
            num_plants = int(width**2*density)
            plant_spacing = width/(num_plants/n_rows[name])



            # Créer une grille de coordonnées pour la catégorie de points actuelle
            grid = []
            for ix in range(n_rows[name]):
                for iy in range(int(num_plants/n_rows[name])):
                    x =  ix * inter_rows[name] + (offset[name] if offset is not None and name in offset else common_offset+2*common_offset*rank)
                    y =  plant_spacing*(iy+0.5)

                    if noise is not None : 
                        (x,y)=(
                            random.uniform(x-noise[name],x+noise[name]),
                            random.uniform(y-noise[name],y+noise[name])
                        )
                    grid.append((x,y,0))
                
            
            rank+=1    #permet de décaler la position des rangs de rank inter-rang (à faire évoluer)
                

            # Stocker la grille de coordonnées dans le dictionnaire
            grids[name] = grid

        #save domain and grid within planter
        domain = ((x_min, y_min), (x_max, y_max))
        coord_grid = grids

        return (domain, coord_grid)
        


    def grid_pattern(self, n_cols= 20, n_rows= 20, cell_size = 0.05, col_pattern = ("wheat","empty","wheat"), noise=None):
        # Créer un dictionnaire pour stocker les grilles de coordonnées de chaque catégorie de points
        grids = {}

        length = n_cols*cell_size
        width = n_rows*cell_size 
        # Calculer les coordonnées min et max du domaine 
        x_min, y_min = 0, 0
        x_max, y_max = length, width

        grids = {name: [] for name in col_pattern}
        
        for name in col_pattern:

            # Créer une grille de coordonnées pour la catégorie de points actuelle
            
            for col_id in range(n_cols):
                grid_col = []
                x =  cell_size/2 + col_id*cell_size
                for row_id in range(n_rows):
                    y =  cell_size/2 + row_id*cell_size
                    if noise is not None : 
                        (x,y)=(
                            random.uniform(x-noise[name],x+noise[name]),
                            random.uniform(y-noise[name],y+noise[name])
                        )
                    grid_col.append((x,y,0))

                modulo =col_id%len(col_pattern)
                name = col_pattern[modulo] 

                # Stocker la grille de coordonnées dans le dictionnaire
                grids[name].extend(grid_col)
            
               

            

        #save domain and grid within planter
        domain = ((x_min, y_min), (x_max, y_max))
        coord_grid = grids

        return (domain, coord_grid)
        

    def __default_preconfigured(
        self, legume_cote={}, inter_rows=0.15, plant_density={1: 250}, xy_plane=None, translate=None, seed=None
    ):
        """Parameters for default generation type

        Note
        ----
        An attribute self.type_domain is created which can be:
            - "create_heterogeneous_canopy": only WheatFspm instance(s) in the simulation, the domain is from AgronomicStand
            - "l-egume": only l-egume instance(s) in the simulation, the domain is from usm configurations
            - "mix": mix of different fspm, we compute global soil domain
            - "input": soil domain == xy_plane

        Parameters
        ----------
        legume_cote : dict, optional
            precise the length of a soil side in l-egume instances. An entry is {wrapper.name : length in cm}, by default {}
        inter_rows : float
            length between two rows in m, by default 0.15
        plant_density : dict, optional
            _description_, by default {1: 250}
        plant_density : dict, optional
            number plants in 1 m^2. An entry is {wrapper.name : nb of plants/m² }, by default {1: 250}
        translate : dict, optional
            Possibility to translate some of the fspm geometric scenes by a 3d vector. An entry is { wrapper.name : (tx, ty, tz) }, by default None
        seed : int, optional
            seed for random and numpy, by default None
        """        
        self.plant_density = plant_density
        self.inter_rows = inter_rows

        for i in self.indexer.wheat_index:
            self.number_of_plants[i] = 50
        
        # runs an iteration of Agronomic stand to get its soil domain
        if self.indexer.wheat_active:
            from alinea.adel.Stand import AgronomicStand
            
            self.type_domain = "create_heterogeneous_canopy"
            self.wheat_positions = [[] for i in range(len(self.indexer.wheat_names))]

            # on vient récupérer le domain de AgronomicStand
            if seed is not None:
                random.seed(seed)
                numpy.random.seed(seed)

            stand = AgronomicStand(
                sowing_density=self.plant_density[1],
                plant_density=self.plant_density[1],
                inter_row=self.inter_rows,
                noise=self.noise_plant_positions,
            )
            _, domain, _, _ = stand.smart_stand(nplants=50, at=self.inter_rows, convunit=1)
            self.domain = domain

        # translate each specy scenes if precised
        if translate is not None:
            self.transformations["translate"] = {}
            for name, vector in translate.items():
                id = self.indexer.global_order.index(name)
                self.transformations["translate"][id] = vector

        # soil domain management
        if xy_plane is None:
            # si recalcul le domain via create_heterogeneous_canopy
            if not self.indexer.legume_active:
                self.type_domain = "create_heterogeneous_canopy"
            elif not self.indexer.wheat_active:
                self.type_domain = "l-egume"
                # convertit domain cm en m
                domains = []
                for name, cote in legume_cote.items():
                    vector = (0, 0, 0)
                    if translate is not None and name in translate :
                        vector = translate[name]
                    domains.append([
                        (0.0 + vector[0], 0.0  + vector[1]),
                        (cote * 0.01  + vector[0], cote * 0.01  + vector[1]),
                    ])
                self.domain = ((min([x[0][0] for x in domains]), min([x[0][1] for x in domains])), 
                            ((max([x[1][0] for x in domains]), max([x[1][1] for x in domains]))))

            else:
                self.type_domain = "mix"
                if self.indexer.legume_active:
                    domains = []
                    for name, cote in legume_cote.items():
                        vector = (0, 0, 0)
                        if translate is not None and name in translate :
                            vector = translate[name]
                        domains.append([
                            (0.0 + vector[0], 0.0  + vector[1]),
                            (cote * 0.01  + vector[0], cote * 0.01  + vector[1]),
                        ])
                    legume_domain = ((min([x[0][0] for x in domains]), min([x[0][1] for x in domains])), 
                                ((max([x[1][0] for x in domains]), max([x[1][1] for x in domains]))))

                # a été calculé au-dessus à  l'appel de create_heterogeneous_canopy
                wheat_domain = self.domain
                if translate is not None:
                    if self.indexer.wheat_names[0] in translate :
                        vector = translate[self.indexer.wheat_names[0]]
                        wheat_domain = (
                            (wheat_domain[0][0] + vector[0], wheat_domain[0][1] + vector[1]),
                            (wheat_domain[1][0] + vector[0], wheat_domain[1][1] + vector[1]),
                        )
                self.domain = (
                    (min(legume_domain[0][0], wheat_domain[0][0]), min(legume_domain[0][1], wheat_domain[0][1])),
                    (max(legume_domain[1][0], wheat_domain[1][0]), max(legume_domain[1][1], wheat_domain[1][1])),
                )
        else:
            self.type_domain = "input"
            self.domain = xy_plane


    def generate_random_other(self, indice_instance=0, seed=None):
        """Compute random plant positions for fspm other than l-egume and wheat

        Parameters
        ----------
        indice_instance : int, optional
            specy ID corresponding to plants to generate, by default 0
        seed : int, optional
            random seed, by default None

        Returns
        -------
        list of tuple
            list of plant positions (x, y, z)
        """        
        if seed is not None:
            s = seed
        else:
            s = 1234
        random.seed(s)
        numpy.random.seed(s)

        # tirage des positions
        # list de 3-tuple des positions
        if self.other_positions[indice_instance] != [] :
            positions = self.other_positions[indice_instance]
        else:
            positions = []
            for i in range(self.number_of_plants[self.indexer.other_index[indice_instance]]):
                positions.append(
                    (numpy.random.uniform(0.0, self.domain[1][0]), numpy.random.uniform(0.0, self.domain[1][0]), 0.0)
                )

        self.other_positions[indice_instance] = positions

        return positions
    
    def generate_random_wheat(
        self, adel_wheat, mtg, indice_wheat_instance=0, stem_name="stem", leaf_name="leaf", seed=None
    ):
        """Compute random plant positions for wheat fspm

        Parameters
        ----------
        adel_wheat : AdelWheat
            adel wheat object containing one wheat geometry
        mtg : openalea.MTG
            MTG of one wheat
        indice_wheat_instance : int, optional
            wheat ID in simulation (indexer.global_index), by default 0
        stem_name : str, optional
            stem tag in the plantgl.Scene, by default "stem"
        leaf_name : str, optional
            leaf tag in the plantgl.Scene, by default "leaf"
        seed : int, optional
            random seed, by default None

        Returns
        -------
        list of tuple
            list of plant positions (x, y, z)
        """        
        var_leaf_inclination = 0.157
        var_leaf_azimut = 1.57
        var_stem_azimut = 0.157

        if seed is not None:
            s = seed
        else:
            s = 1234
        random.seed(s)
        numpy.random.seed(s)

        initial_scene = adel_wheat.scene(mtg)

        # tirage des positions
        # list de 3-tuple des positions
        if self.wheat_positions[indice_wheat_instance] != [] and self.save_wheat_positions:
            positions = self.wheat_positions[indice_wheat_instance]
        else:
            positions = []
            for i in range(self.number_of_plants[self.indexer.wheat_index[indice_wheat_instance]]):
                positions.append(
                    (numpy.random.uniform(0.0, self.domain[1][0]), numpy.random.uniform(0.0, self.domain[1][0]), 0.0)
                )

        self.wheat_positions[indice_wheat_instance] = positions

        generated_scene = self.__generate_wheat_from_positions(
            initial_scene, mtg, positions, var_leaf_inclination, var_leaf_azimut, var_stem_azimut, stem_name, leaf_name
        )

        return generated_scene

    def generate_row_other(self, indice_instance=0, seed=None):
        """Compute row plant positions for fspm other than l-egume and wheat

        Parameters
        ----------
        indice_instance : int, optional
            specy ID corresponding to plants to generate, by default 0
        seed : int, optional
            random seed, by default None

        Returns
        -------
        list of tuple
            list of plant positions (x, y, z)
        """        

        if seed is not None:
            s = seed
        else:
            s = 1234
        random.seed(s)
        numpy.random.seed(s)

        if self.other_positions[indice_instance] != []:
            positions = self.other_positions[indice_instance]
        else:
            positions = []

            inter_plants = (#space between plants within the row
                2 * self.domain[1][1] / self.number_of_plants[self.indexer.wheat_index[indice_instance]]
                #2 * length of the row / number of plants in this fspm instance, as indicated by the indexer
            )
            nrows = 2
            self.total_n_rows

            
            rows_y = [self.inter_rows * 1.5, ((self.total_n_rows / nrows) + 1.5) * self.inter_rows]
            #rows y coordinates = [1st row is placed at a distance of 1.5 inter row from the edge of the domain, 2nd row is placed X inter-rows further, X being the ratio ]
            for y in rows_y:
                for ix in range(int(self.number_of_plants[self.indexer.wheat_index[indice_instance]] / nrows)):
                    x = inter_plants * (0.5 + ix)
                    p = (
                        random.uniform(x - self.noise_plant_positions, x + self.noise_plant_positions),
                        random.uniform(y - self.noise_plant_positions, y + self.noise_plant_positions),
                        0.0,
                    )
                    positions.append(p)

        self.other_positions[indice_instance] = positions

        return positions

    def generate_row_wheat(
        self, adel_wheat, mtg, indice_wheat_instance=0, stem_name="stem", leaf_name="leaf", seed=None
    ):
        """Compute row plant positions for wheat fspm

        Parameters
        ----------
        adel_wheat : AdelWheat
            adel wheat object containing one wheat geometry
        mtg : openalea.MTG
            MTG of one wheat
        indice_wheat_instance : int, optional
            wheat ID in simulation (indexer.global_index), by default 0
        stem_name : str, optional
            stem tag in the plantgl.Scene, by default "stem"
        leaf_name : str, optional
            leaf tag in the plantgl.Scene, by default "leaf"
        seed : int, optional
            random seed, by default None

        Returns
        -------
        list of tuple
            list of plant positions (x, y, z)
        """        
        var_leaf_inclination = 0.157
        var_leaf_azimut = 1.57
        var_stem_azimut = 0.157

        if seed is not None:
            s = seed
        else:
            s = 1234
        random.seed(s)
        numpy.random.seed(s)

        initial_scene = adel_wheat.scene(mtg)

        if self.wheat_positions[indice_wheat_instance] != [] and self.save_wheat_positions:
            positions = self.wheat_positions[indice_wheat_instance]
        else:
            positions = []

            inter_plants = ( #space between plants within the row
                2 * self.domain[1][1] / self.number_of_plants[self.indexer.wheat_index[indice_wheat_instance]]
                #2 * length of the row / number of plants in this fspm instance, as indicated by the indexer
            )
            nrows = 2
            self.total_n_rows

            # first row on left 1/2 interrow, then 1 out of 2 row is wheat
            rows_y = [self.inter_rows * 1.5, ((self.total_n_rows / nrows) + 1.5) * self.inter_rows]
            # y coordinates of the rows = 
            for y in rows_y:
                for ix in range(int(self.number_of_plants[self.indexer.wheat_index[indice_wheat_instance]] / nrows)):
                    x = inter_plants * (0.5 + ix)
                    p = (
                        random.uniform(x - self.noise_plant_positions, x + self.noise_plant_positions),
                        random.uniform(y - self.noise_plant_positions, y + self.noise_plant_positions),
                        0.0,
                    )
                    positions.append(p)

        self.wheat_positions[indice_wheat_instance] = positions

        generated_scene = self.__generate_wheat_from_positions(
            initial_scene,
            mtg,
            positions,
            var_leaf_inclination,
            var_leaf_azimut,
            var_stem_azimut,
            stem_name=stem_name,
            leaf_name=leaf_name,
        )

        return generated_scene

    def generate_forced_wheat(
            self, adel_wheat, mtg, indice_wheat_instance=0, stem_name="stem", leaf_name="leaf", seed=None
    ):
        
        """Generates wheat canopy based on position grid

        Parameters
        ----------
        adel_wheat : AdelWheat
            adel wheat object containing one wheat geometry
        mtg : openalea.MTG
            MTG of one wheat
        indice_wheat_instance : int, optional
            wheat ID in simulation (indexer.global_index), by default 0
        stem_name : str, optional
            stem tag in the plantgl.Scene, by default "stem"
        leaf_name : str, optional
            leaf tag in the plantgl.Scene, by default "leaf"
        seed : int, optional
            random seed, by default None

        Returns
        -------
        a generated scene 
        """        
        var_leaf_inclination = 0.157
        var_leaf_azimut = 1.57
        var_stem_azimut = 0.157
        if seed is not None:
            s = seed
        else:
            s = 1234
        random.seed(s)
        numpy.random.seed(s)

        self.save_wheat_positions=True

        initial_scene = adel_wheat.scene(mtg)

        positions = self.wheat_positions[indice_wheat_instance]

        generated_scene = self.__generate_wheat_from_positions(
            initial_scene,
            mtg,
            positions,
            var_leaf_inclination,
            var_leaf_azimut,
            var_stem_azimut,
            stem_name=stem_name,
            leaf_name=leaf_name,
        )

        return generated_scene






    def create_heterogeneous_canopy(
        self,
        geometrical_model,
        mtg=None,
        var_leaf_inclination=0.157,
        var_leaf_azimut=1.57,
        var_stem_azimut=0.157,
        stem_name="stem",
        leaf_name="leaf",
        indice_wheat_instance=0,
        seed=None,
    ):
        """Generate wheat positions in default mode

        Parameters
        ----------
        geometrical_model : AdelWheat
            adel wheat object the wheat
        mtg : openalea.MTG, optional
            MTG containing plantgl.Scene of one wheat, by default None
        var_leaf_inclination : float, optional
            variability for leaf inclination (rad), by default 0.157
        var_leaf_azimut : float, optional
            variability for leaf azimut (rad), by default 1.57
        var_stem_azimut : float, optional
            variability for stem azimut (rad), by default 0.157
        stem_name : str, optional
            stem tag in the plantgl.Scene, by default "stem"
        leaf_name : str, optional
            leaf tag in the plantgl.Scene, by default "leaf"
        indice_wheat_instance : int, optional
            wheat specy ID in simulation, by default 0
        seed : int, optional
            random seed, by default None

        Returns
        -------
        list of tuple
            list of plant positions (x, y, z)
        """        

        from alinea.adel.Stand import AgronomicStand
        
        if seed is not None:
            random.seed(seed)
            numpy.random.seed(seed)

        # Planter
        stand = AgronomicStand(
            sowing_density=self.plant_density[1],
            plant_density=self.plant_density[1],
            inter_row=self.inter_rows,
            noise=self.noise_plant_positions,
        )
        _, domain, positions, _ = stand.smart_stand(
            nplants=self.number_of_plants[self.indexer.wheat_index[indice_wheat_instance]],
            at=self.inter_rows,
            convunit=1,
        )
        self.wheat_positions[indice_wheat_instance] = positions

        random.seed(1234)

        generated_scene = self.__generate_wheat_from_positions(
            geometrical_model,
            mtg,
            positions,
            var_leaf_inclination,
            var_leaf_azimut,
            var_stem_azimut,
            stem_name=stem_name,
            leaf_name=leaf_name,
        )

        if self.type_domain == "create_heterogeneous_canopy":
            self.domain = domain

        return generated_scene

    def __generate_wheat_from_positions(
        self,
        geometrical_model,
        mtg=None,
        positions=[(0.0, 0.0, 0.0)],
        var_leaf_inclination=0.157,
        var_leaf_azimut=1.57,
        var_stem_azimut=0.157,
        stem_name="stem",
        leaf_name="leaf",
    ):
        """Generate complete wheats from their plant positions

        Parameters
        ----------
        geometrical_model : AdelWheat
            adel wheat object the wheat
        mtg : openalea.MTG, optional
            MTG containing plantgl.Scene of one wheat, by default None
        positions : list, optional
            list of plant positions as (x, y, z) in m, by default [(0.0, 0.0, 0.0)]
        var_leaf_inclination : float, optional
            variability for leaf inclination (rad), by default 0.157
        var_leaf_azimut : float, optional
            variability for leaf azimut (rad), by default 1.57
        var_stem_azimut : float, optional
            variability for stem azimut (rad), by default 0.157
        stem_name : str, optional
            stem tag in the plantgl.Scene, by default "stem"
        leaf_name : str, optional
            leaf tag in the plantgl.Scene, by default "leaf"

        Returns
        -------
        plantgl.Scene
            final scene of the generated wheats 
        """        
        # Load scene
        if not isinstance(geometrical_model, plantgl.Scene):
            initial_scene = geometrical_model.scene(mtg)
        else:
            initial_scene = geometrical_model

        alea_canopy = pandas.DataFrame()

        # Built alea table if does not exist yet
        if alea_canopy.empty and mtg is not None:
            elements_vid_list = []
            for mtg_plant_vid in mtg.components_iter(mtg.root):
                for mtg_axis_vid in mtg.components_iter(mtg_plant_vid):
                    for mtg_metamer_vid in mtg.components_iter(mtg_axis_vid):
                        for mtg_organ_vid in mtg.components_iter(mtg_metamer_vid):
                            for mtg_element_vid in mtg.components_iter(mtg_organ_vid):
                                if mtg.label(mtg_element_vid) == leaf_name:
                                    elements_vid_list.append(mtg_element_vid)

            elements_vid_df = pandas.DataFrame({"vid": elements_vid_list, "tmp": 1})
            positions_df = pandas.DataFrame(
                {"pos": range(len(positions)), "tmp": 1, "azimut_leaf": 0, "inclination_leaf": 0}
            )
            alea = pandas.merge(elements_vid_df, positions_df, on=["tmp"])
            alea = alea.drop("tmp", axis=1)
            for vid in elements_vid_list:
                numpy.random.seed(vid)
                alea.loc[alea["vid"] == vid, "azimut_leaf"] = numpy.random.uniform(
                    -var_leaf_azimut, var_leaf_azimut, size=len(positions)
                )
                alea.loc[alea["vid"] == vid, "inclination_leaf"] = numpy.random.uniform(
                    -var_leaf_inclination, var_leaf_inclination, size=len(positions)
                )
            alea_canopy = alea

        # Duplication and heterogeneity
        duplicated_scene = plantgl.Scene()
        position_number = 0
        for pos in positions:
            azimut_stem = random.uniform(-var_stem_azimut, var_stem_azimut)
            for shp in initial_scene:
                if mtg.label(shp.id) == stem_name:
                    rotated_geometry = plantgl.EulerRotated(azimut_stem, 0, 0, shp.geometry)
                    translated_geometry = plantgl.Translated(plantgl.Vector3(pos), rotated_geometry)
                    new_shape = plantgl.Shape(translated_geometry, appearance=shp.appearance, id=shp.id)
                    duplicated_scene += new_shape
                elif mtg.label(shp.id) == leaf_name:
                    # Add shp.id in alea_canopy if not in yet:
                    if shp.id not in list(alea_canopy["vid"]):
                        new_vid_df = pandas.DataFrame({"vid": shp.id, "pos": range(len(positions))})
                        numpy.random.seed(shp.id)
                        new_vid_df["azimut_leaf"] = numpy.random.uniform(
                            -var_leaf_azimut, var_leaf_azimut, size=len(positions)
                        )
                        new_vid_df["inclination_leaf"] = numpy.random.uniform(
                            -var_leaf_inclination, var_leaf_inclination, size=len(positions)
                        )
                        alea_canopy = alea_canopy.copy().append(new_vid_df, sort=False)
                    # Translation to origin
                    anchor_point = mtg.get_vertex_property(shp.id)["anchor_point"]
                    trans_to_origin = plantgl.Translated(-anchor_point, shp.geometry)
                    # Rotation variability
                    azimut = alea_canopy.loc[
                        (alea_canopy.pos == position_number) & (alea_canopy.vid == shp.id), "azimut_leaf"
                    ].values[0]
                    inclination = alea_canopy.loc[
                        (alea_canopy.pos == position_number) & (alea_canopy.vid == shp.id), "inclination_leaf"
                    ].values[0]
                    rotated_geometry = plantgl.EulerRotated(azimut, inclination, 0, trans_to_origin)
                    # Restore leaf base at initial anchor point
                    translated_geometry = plantgl.Translated(anchor_point, rotated_geometry)
                    # Translate leaf to new plant position
                    translated_geometry = plantgl.Translated(pos, translated_geometry)
                    new_shape = plantgl.Shape(translated_geometry, appearance=shp.appearance, id=shp.id)
                    duplicated_scene += new_shape

            position_number += 1

        return duplicated_scene
