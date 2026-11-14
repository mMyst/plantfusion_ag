import os
import zipfile
import tkinter as tk
from tkinter import filedialog

def select_directory() -> str:
    """Opens a native folder selection dialog and returns the selected path."""
    root = tk.Tk()
    root.withdraw()  # Hide the default Tk window
    root.attributes("-topmost", True)  # Bring the dialog to the front

    folder_selected = filedialog.askdirectory(
        title="Select Folder to Crawl for Zip Files"
    )
    root.destroy()
    return folder_selected

def extract_and_remove_zips(root_dir: str):
    """
    Recursively finds all .zip files within root_dir,
    extracts them into their containing directory, and deletes the original .zip.
    """
    zip_files = []
    for dirpath, _, filenames in os.walk(root_dir):
        for filename in filenames:
            if filename.lower().endswith(".zip"):
                zip_files.append(os.path.join(dirpath, filename))

    if not zip_files:
        print("No zip files found.")
        return

    print(f"Found {len(zip_files)} zip file(s). Processing...")

    for zip_path in zip_files:
        extract_dir = os.path.dirname(zip_path)
        try:
            with zipfile.ZipFile(zip_path, "r") as archive:
                archive.extractall(extract_dir)
            
            # Delete only after successful extraction
            os.remove(zip_path)
            print(f"Extracted and deleted: {zip_path}")
            
        except zipfile.BadZipFile:
            print(f"Skipped (corrupt/invalid archive): {zip_path}")
        except PermissionError:
            print(f"Permission denied: {zip_path}")
        except Exception as e:
            print(f"Failed to process {zip_path}: {e}")

if __name__ == "__main__":
    target_directory = select_directory()

    if target_directory:
        print(f"Selected: {target_directory}")
        extract_and_remove_zips(target_directory)
    else:
        print("Operation cancelled: No folder was selected.")