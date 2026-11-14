import os
import numpy as np
import openalea.plantgl.all as pgl

def get_top_triangles_height(scene, percentage=10, height_type='barycenter'):
    """
    Calculates the average height of the top X% highest triangles in a PlantGL scene.
    
    :param scene: The PlantGL Scene object containing the 3D shapes.
    :param percentage: The percentage of the highest triangles to consider (default 10).
    :param height_type: 'barycenter' (mean Z of the 3 vertices) or 'max' (highest vertex Z).
    :return: The average height of the top X% triangles.
    """
    triangle_heights = []
    
    # Initialize a Tesselator to convert any non-triangle shapes into triangles
    # tesselator = pgl.Tesselator()
    
    # Iterate through all shapes in the PlantGL scene
    for shape in scene:
        geom = shape.geometry
        
        # # If the geometry is not already a TriangleSet, try to triangulate it
        # if not isinstance(geom, pgl.TriangleSet):
        #     try:
        #         geom.apply(tesselator)
        #         geom = tesselator.result
        #     except Exception as e:
        #         # Skip if the shape cannot be triangulated
        #         print(f"Warning: Could not triangulate shape {shape.id}: {e}")
        #         continue
                
        # Process the TriangleSet
        if isinstance(geom, pgl.TriangleSet):
            points = geom.pointList
            indices = geom.indexList
            
            # Iterate through each triangle
            for idx in indices:
                # Extract the Z coordinate (upward axis) for each vertex
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
        print("Warning: No triangles found in the scene.")
        return 0.0
        
    # Find the threshold for the top X%
    threshold = np.percentile(heights_array, 100 - percentage)
    
    # Filter and calculate the mean
    top_heights = heights_array[heights_array >= threshold]
    
    return float(np.mean(top_heights))


# ==========================================
# Testing the function with a .bgeom file
# ==========================================
if __name__ == "__main__":
    
    # Replace this with the actual path to your .bgeom file
    bgeom_file_path = r"C:\Users\agrumel\Code\Python_Ecophy\plantfusion_ag\outputs\WheatLuz_LUBBAC_emergence_adel_autoreg\light\plantgl\scene_light_plantgl_10333513.bgeom"
    
    if not os.path.exists(bgeom_file_path):
        print(f"Error: The file '{bgeom_file_path}' does not exist.")
        print("Please update 'bgeom_file_path' with a valid path to test the script.")
    else:
        print(f"Loading scene from '{bgeom_file_path}'...")
        
        # Load the bgeom file into a PlantGL Scene
        my_scene = pgl.Scene(bgeom_file_path)
        
        # Test 1: Calculate using the barycenter (center of the triangle)
        top_10_barycenter = get_top_triangles_height(my_scene, percentage=10, height_type='barycenter')
        print(f"Average height of the top 10% triangles (barycenter): {top_10_barycenter:.2f} m")
        
        # Test 2: Calculate using the maximum vertex (absolute highest point of the triangle)
        top_10_max = get_top_triangles_height(my_scene, percentage=10, height_type='max')
        print(f"Average height of the top 10% triangles (highest vertex): {top_10_max:.2f} m")
        
        # Test 3: Change the percentage (e.g., top 5%)
        top_5_barycenter = get_top_triangles_height(my_scene, percentage=5, height_type='barycenter')
        print(f"Average height of the top 5% triangles (barycenter): {top_5_barycenter:.2f} m")