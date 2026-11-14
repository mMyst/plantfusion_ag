import os
import warnings

import numpy as np
from PIL import Image, ImageFilter, ImageOps
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import openalea.plantgl.all as pgl


def save_df_to_csv(df, outputs_filepath, precision):
    """Save pandas dataframes to csv

    Parameters
    ----------
    df : pandas.DataFrame
        a pandas dataframe to be saved
    outputs_filepath : str
        the name of the CSV file to be saved
    precision : int
        number of decimals in CSV file
    """    

    try:
        df.to_csv(outputs_filepath, na_rep="NA", index=False, float_format="%.{}f".format(precision))
    except IOError as err:
        path, filename = os.path.split(outputs_filepath)
        filename = os.path.splitext(filename)[0]
        newfilename = "ACTUAL_{}.csv".format(filename)
        newpath = os.path.join(path, newfilename)
        df.to_csv(newpath, na_rep="NA", index=False, float_format="%.{}f".format(precision))
        warnings.warn("[{}] {}".format(err.errno, err.strerror))
        warnings.warn("File will be saved at {}".format(newpath))


def create_child_folder(parentfolderpath, childfolderpath):
    """Create a child folder from parent folder path

    Parameters
    ----------
    parentfolderpath : string
        parent folder path
    childfolderpath : string
        name of the child folder
    """    
    dirName = os.path.join(os.path.normpath(parentfolderpath), os.path.normpath(childfolderpath))
    try:
        # Create target Directory
        os.mkdir(dirName)
        print("Directory ", dirName, " Created ")
    except FileExistsError:
        print("Directory ", dirName, " already exists")

def convert_colors(image, n_colors):
    # Convert image to NumPy array
    image_np = np.array(image)

    # Reduce the number of image dimensions
    image_np_reshape = image_np.reshape(-1, 3)

    # Apply KMeans algorithm
    kmeans = KMeans(n_clusters=n_colors)
    labels = kmeans.fit_predict(image_np_reshape)

    # Create a new image from the dominant colors
    image_kmeans = kmeans.cluster_centers_[labels].reshape(image_np.shape)

    # Convert NumPy array to PIL image
    image_resultat = Image.fromarray(image_kmeans.astype('uint8'))

    return image_resultat, image_kmeans

def convert_nb(image, n_nuances):
    # Convert image to grayscale
    image_grise = image.convert('L')

    # Convert grayscale image to NumPy array
    image_np = np.array(image_grise)

    # Calculate quantization step
    pas = 256 / n_nuances

    # Quantize the image using the calculated step
    image_quantifiee = np.digitize(image_np, np.arange(0, 256, pas)) * pas

    # Convert NumPy array to PIL image
    image_resultat = Image.fromarray(image_quantifiee.astype('uint8'))

    return image_resultat, image_quantifiee

def convert_portrait(image, threshold=50):
    """Trace les contours d'un visage pour un mode portrait."""
    # Convertir en niveaux de gris
    image_grise = image.convert('L')
    
    # Appliquer le filtre FIND_EDGES pour détecter les contours
    edges = image_grise.filter(ImageFilter.FIND_EDGES)
    
    # Inverser les couleurs pour obtenir des traits noirs sur un fond blanc
    edges_inverted = ImageOps.invert(edges)
    
    # Binariser l'image (Nettoyage pour n'avoir que du noir pur et du blanc pur)
    image_np = np.array(edges_inverted)
    image_binaire = np.where(image_np > threshold, 255, 0).astype('uint8')
    
    # Reconvertir le tableau NumPy en image PIL
    image_resultat = Image.fromarray(image_binaire)
    
    return image_resultat, image_binaire


def coord_from_image(image, n_colors, n_side, convert_mode=None,scale=10**-2):
    """
    Cuts an image into a grid of n squares, changes the color of each square to the most abundant color in that square,
    and returns the modified image as well as the coordinates of the barycenters of each square.

    Args:
        image (PIL.Image): The image to be cut.
        n (int): The number of squares per side of the grid.

    Returns:
        tuple: A tuple containing the modified image and a dictionary of barycenter coordinates for each square.
    """

    n = n_side

    if convert_mode is not None:
        if convert_mode == 'nb':
            image, _ = convert_nb(image, n_colors) #
        elif convert_mode == 'colors':
            image, _ = convert_colors(image, n_colors) #
        elif convert_mode == 'portrait':
            # Appel de la nouvelle fonction portrait
            image, _ = convert_portrait(image) 
        else:
            print("Unknown conversion mode.") #

    # Convert the image to a NumPy array and flip it 
    array = np.array(image)
    array = np.fliplr(array)
    # Calculate the size of each square
    height, width = array.shape
    size = max(height, width) // n
    n_long = int(height / size)
    n_larg = int(width / size)
    # Create an empty image of the same size as the original image
    image_modifiee = Image.new("RGB", image.size)
    # Create a dictionary to store the coordinates of the barycenters of each square
    coords = {color: [] for color in np.unique(array)}
    # Loop through each square of the grid
    for i in range(n_long):
        for j in range(n_larg):
            # Calculate the coordinates of the square in the original image
            x1 = j * size
            y1 = i * size
            x2 = min(width, (j + 1) * size)
            y2 = min(height, (i + 1) * size)
            
            # Extract the square from the original image
            carre = array[y1:y2, x1:x2]
            
            # --- NOUVELLE LOGIQUE POUR LA COULEUR ---
            if convert_mode == 'portrait':
                # Compter le nombre de pixels foncés (noirs) dans le carré
                pixels_noirs = np.sum(carre < 128)
                # Si plus de 5% du carré contient un trait, on force le carré en noir (0)
                seuil = (carre.shape[0] * carre.shape[1]) * 0.05
                color = 0 if pixels_noirs > seuil else 255
            else:
                # Calculate the most abundant color in the square (Mode original)
                color = np.bincount(carre.ravel(), minlength=256).argmax()
            # ----------------------------------------

            # Change the color of the square to the selected color
            image_modifiee.paste(int(color), (x1, y1, x2, y2))
            
            # Calculate the coordinates of the barycenter of the square
            x = (x1 + (x2 - x1) // 2)*scale
            y = (y1 + (y2 - y1) // 2)*scale
            
            # Add the barycenter coordinates to the dictionary
            coords[color].append((x, y, 0))
            
    return image_modifiee, coords

def get_top_triangles_height(scene, percentage=10, height_type='barycenter'):
    """
    Calculates the average height of the top X% highest triangles in a PlantGL scene.
    
    :param scene: The PlantGL Scene object containing the 3D shapes.
    :param percentage: The percentage of the highest triangles to consider (default 10).
    :param height_type: 'barycenter' (mean Z of the 3 vertices) or 'max' (highest vertex Z).
    :return: The average height of the top X% triangles.
    """
    triangle_heights = []
    
    # Iterate through all shapes in the PlantGL scene
    for shape in scene:
        geom = shape.geometry
        
        # Check if the geometry is made of triangles
        # Note: If your scene uses other shapes (cylinders, spheres), 
        # they must be triangulated first (e.g., using pgl.Tesselator)
        if isinstance(geom, pgl.TriangleSet):
            points = geom.pointList
            indices = geom.indexList
            
            # Iterate through each triangle
            for idx in indices:
                # idx contains the 3 indices of the vertices in the pointList
                # We extract the Z coordinate (upward axis) for each vertex
                z1 = points[idx[0]].z
                z2 = points[idx[1]].z
                z3 = points[idx[2]].z
                
                # Define the height of the current triangle
                if height_type == 'barycenter':
                    tri_height = (z1 + z2 + z3) / 3.0
                else:
                    tri_height = max(z1, z2, z3)
                    
                triangle_heights.append(tri_height)
                
    # Convert the list to a NumPy array for efficient calculation
    heights_array = np.array(triangle_heights)
    
    # Safety check if no triangles were found
    if heights_array.size == 0:
        return 0.0
        
    # Find the threshold for the top X%
    threshold = np.percentile(heights_array, 100 - percentage)
    
    # Filter and calculate the mean
    top_heights = heights_array[heights_array >= threshold]
    
    return float(np.mean(top_heights))

def get_top_height(scene, percentage=10):
    """
    Calculates the average height of the top X% highest shapes in a PlantGL scene.
    Keeps things simple by using PlantGL's native BoundingBox.
    """
    heights = []
    
    # Iterate through every object (leaf, stem, etc.) in the scene
    for shape in scene:
        # PlantGL automatically handles all transformations (Translated, etc.) here
        bbox = pgl.BoundingBox(shape)
        
        # .upper returns a Vector3 (x, y, z) of the absolute highest point of this shape
        max_z = bbox.upperRightCorner.z
        
        heights.append(max_z)
        
    heights_array = np.array(heights)
    
    # Safety check if the scene is empty
    if heights_array.size == 0:
        return 0.0
        
    # Find the threshold for the top X%
    threshold = np.percentile(heights_array, 100 - percentage)
    
    # Filter and calculate the mean
    top_heights = heights_array[heights_array >= threshold]
    
    return float(np.mean(top_heights))