import pandas as pd
import shutil

def remove_folder(folder_path):
    try:
        shutil.rmtree(folder_path)
        print(f"Folder '{folder_path}' removed successfully.")
    except FileNotFoundError:
        print(f"Folder '{folder_path}' not found.")
    except PermissionError:
        print(f"Permission denied to delete folder '{folder_path}'.")


input_file = "../Database_campioni2.csv"
output_file = "../Database_campioni_corretto.csv"

try:
    df = pd.read_csv(input_file, index_col=0)
except FileNotFoundError:
    print(f"Error: {input_file} not found.")
    exit()

# Step 2: Modify all the values in the column "seed" with incremental integer numbers
df_no_duplicates = df.drop_duplicates()
duplicated_rows = df.loc[df.duplicated()]
#for i in duplicated_rows.index:
#    remove_folder(f"../samples_2/{i}")
duplicated_rows['Seed'] = range(-len(duplicated_rows), 0)
corrected_df = pd.concat([df_no_duplicates, duplicated_rows])
# Step 3: Write the dataframe into a new .csv file
corrected_df.to_csv(output_file)

print(f"Modified dataframe saved to {output_file}.")
