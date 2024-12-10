import pandas as pd
import matplotlib.pyplot as plt

campioni = pd.read_csv("Database_campioni2.csv", index_col=0)
t = campioni[campioni["Size"] == 0.1]
#Selecting train samples
t = t[t["Real/Simul"] == "Simul"]
Train = t[t["Aureus_S_N_strains"] == 5]
Train_abundances = Train[["Aureus_S", "Aureus_R", "Coli", "Pseudomonas", "Arcobacter"]]

#plot single species abundance in the train set
sums = Train_abundances.sum()
plt.figure(figsize=(10, 6))
sums.plot(kind='bar')
plt.title('Somma delle abbondanze per ciascuna specie nel set di addestramento')
plt.xlabel('Specie')
plt.ylabel('Somma delle abbondanze')
plt.xticks(rotation=45)
plt.show()

#plot the number of samples in the train set that contain each species
Train_present = Train_abundances[Train_abundances != 0].count()
plt.figure(figsize=(10, 6))
Train_present.plot(kind='bar')
plt.title('Numero di campioni nel set di addestramento che contengono ciascuna specie')
plt.xlabel('Specie')
plt.ylabel('Numero di campioni')
plt.xticks(rotation=45)
plt.show()

