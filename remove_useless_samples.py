import pandas as pd
import os
import shutil
#Lunghezza media delle reads simulate e reali, comunque un fattore di normalizzazione che dovrebbe valere per entrambi
#potrebbe essere diverso tra simulate e reali
average_read_length = 270
images_dir = "M_learning_images/"
Database_name = "Database_campioni2.csv"
campioni = pd.read_csv(Database_name, index_col=0)

subsample_counts = 1
USE_BIAS = 2
WHICH = "_"

#Deciding which samples to peack from the database and defining the train and test sets
t = campioni[campioni["Size"] == 0.01]
t = t.index
t = t.tolist()
path_t = ["/home/PERSONALE/adriano.zaghi2/AMR/samples_2/" + i for i in t]

i=0
for dir_path in path_t:
    if os.path.exists(dir_path) and os.path.isdir(dir_path):
        #shutil.rmtree(dir_path)
        #print(f"Removed directory: {dir_path}")
        i += 1
    else:
        print(f"Directory does not exist or is not a directory: {dir_path}")

print(i)
cleaned_database = campioni.drop(t)
cleaned_database.to_csv("Clean_"+Database_name)