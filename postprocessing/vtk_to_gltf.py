"""
Script to convert VTK files to GLTF format using LightVegeManager.
Usage: python vtk_to_gltf.py <input_vtk> <output_gltf>
"""

import sys
import os
from lightvegemanager.GLTF import VTKtoGLTF

def main():
    if len(sys.argv) != 3:
        print("Usage: python vtk_to_gltf.py <input_vtk_file> <output_gltf_file>")
        return

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' does not exist.")
        return

    print(f"Converting '{input_file}' to '{output_file}'...")
    
    try:
        VTKtoGLTF(input_file, output_file)
        print("Conversion completed successfully.")
    except Exception as e:
        print(f"An error occurred during conversion: {e}")

if __name__ == "__main__":
    main()