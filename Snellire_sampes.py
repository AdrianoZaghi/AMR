import os
import shutil

# Define the base directory where all structured folders are located
BASE_DIR = "samples_4"

# Walk through all subdirectories inside the base directory
for root, dirs, files in os.walk(BASE_DIR):
    item_path = os.path.join(root, "camisim_out")
    #check if the directory actually exists
    if os.path.isdir(item_path):
        shutil.rmtree(item_path)

print("Cleanup complete!")
