import os
import warnings

import numpy as np
from PIL import Image
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt


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
            image, _ = convert_nb(image, n_colors)
        elif convert_mode == 'colors':
            image, _ = convert_colors(image, n_colors)
        else:
            print("Unknown conversion mode.")

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
            # Calculate the most abundant color in the square
            color = np.bincount(carre.ravel(), minlength=256).argmax()
            # Change the color of the square to the most abundant color
            image_modifiee.paste(int(color), (x1, y1, x2, y2))
            # Calculate the coordinates of the barycenter of the square
            x = (x1 + (x2 - x1) // 2)*scale
            y = (y1 + (y2 - y1) // 2)*scale
            # Add the barycenter coordinates to the dictionary
            coords[color].append((x, y, 0))
    return image_modifiee, coords