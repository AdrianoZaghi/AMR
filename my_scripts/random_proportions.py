import random
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def genera_proporzioni(Batteri, lista = ["Aureus_S","Aureus_R","Coli","Pseudomonas","Arcobacter"]):
    if len([x for x in Batteri if x not in lista]) > 0:
        raise ValueError("Batteri non validi. Devono essere tra:" + str(lista) + ".")
    quanti_te_ne_servono = len(Batteri)
    is_ok = False
    while is_ok == False:
        #numbers = sorted([random.randint(0,1000) for x in range(0, quanti_te_ne_servono-1)])
        #numbers.append(1000)
        numbers = sorted([random.randint(0,1000) for x in range(0, quanti_te_ne_servono)])
        max = numbers[-1]
        if max == 0:
            continue
        numbers.insert(0,0)
        numbers = [100*x/max for x in numbers]
        intervals = [(numbers[i+1] - numbers[i]) for i in range(0,len(numbers)-1)]
        if sum(intervals) == 100:
                is_ok = True
    to_return = []
    for i in lista:
        if i in Batteri:
            to_return.append(intervals.pop(0))
        else:
            to_return.append(0)
    return to_return

#Aureus_S,Aureus_R,Coli,Pseudomonas,Arcobacter
#Arcobacter è quello suscettibile

def genera_df(Quanti_campioni, Pronome, Batteri, Size, N_strains, Lower_Seed = 1):
    data = []
    for i in range(0,Quanti_campioni):
        Nome = Pronome + str(i)
        proporzioni = genera_proporzioni(Batteri)
        sample = [Nome] + proporzioni + [Size, N_strains, "De_Caesar_simul", "Campione_simulato", Lower_Seed + 1 + i]
        data.append(sample)
    df = pd.DataFrame(data, columns=["Nome","Aureus_S","Aureus_R","Coli","Pseudomonas","Arcobacter","Size","N strains","Real/Simul","Original sample","Seed"])
    return df

def genera_df_2(Quanti_campioni, Pronome, Batteri, Size, Strain_counts, Lower_Seed = 1):
    data = []
    if "Aureus_S" not in Batteri:
        Strain_counts.insert(0,0)
    if "Aureus_R" not in Batteri:
        Strain_counts.insert(1,0)
    if "Coli" not in Batteri:
        Strain_counts.insert(2,0)
    if "Pseudomonas" not in Batteri:
        Strain_counts.insert(3,0)
    if "Arcobacter" not in Batteri:
        Strain_counts.insert(4,0)
    for i in range(0,Quanti_campioni):
        Nome = Pronome + str(i)
        proporzioni = genera_proporzioni(Batteri)
        sample = [Nome] + proporzioni + [Size, "De_Caesar_simul", "Campione_simulato", Lower_Seed + 1 + i] + Strain_counts
        print(sample)
        data.append(sample)
    df = pd.DataFrame(data, columns=["Nome","Aureus_S","Aureus_R","Coli","Pseudomonas","Arcobacter","Size","Real/Simul","Original sample","Seed","SeedAureus_S_N_strains","Aureus_R_N_strains","Coli_N_strains","Pseudomonas_N_strains","Arcobacter_N_strains"])
    return df

def genera_df_3(Quanti_campioni, Pronome, Size, Strain_counts, Genome_source : ["Simul", "De_Caesar_simul"], Lower_Seed = 1):
    """Generate "Quanti_campioni" samples with name "Pronome" + i, with size "Size" and strain counts "Strain_counts" and seed "Seed" specified in the arguments.
    The presence of a strain is chosen randomly, first selecting how many of them, and then selecting them from our genome pool.
    The number of symulated strains is fixed for each of the five strains, and is specified in the "Strain_counts" list.
    Specifying a number of strains different from 0 for a strain not present in the sample, will not be a problem: if the bacteria abundance is selected to be 0, nothing of it will go in to the sample.
    """
    data = []
    possible_bacteria = ["Aureus_S","Aureus_R","Coli","Pseudomonas","Arcobacter"]
    if len(Strain_counts) != 5:
        raise ValueError("Strain_counts deve avere lunghezza 5")
    for i in range(0,Quanti_campioni):
        quanti_batteri = random.randint(1,5)
        quali_batteri = random.sample(possible_bacteria, quanti_batteri)
        pop = genera_proporzioni(quali_batteri)
        Nome = Pronome + str(i)
        provenienza = Genome_source
        sample = [Nome] + pop + [Size, provenienza, "Campione_simulato", Lower_Seed + i + 1] + Strain_counts
        data.append(sample)
    df = pd.DataFrame(data, columns=["Nome"] + possible_bacteria +["Size","Real/Simul","Original sample","Seed","Aureus_S_N_strains","Aureus_R_N_strains","Coli_N_strains","Pseudomonas_N_strains","Arcobacter_N_strains"])
    df.set_index("Nome", inplace = True)
    return df

if False:
    #controllando la distribuzione delle proporzioni generate dal mio modello
    data = []
    for i in range(0, 10000000):
        data.append(genera_proporzioni(["Aureus_S","Aureus_R","Coli","Pseudomonas","Arcobacter"]))
    df = pd.DataFrame(data, columns=["Aureus_S","Aureus_R","Coli","Pseudomonas","Arcobacter"])
    colors = ['red', 'blue', 'green', 'purple', 'orange']
    # Plot a histogram for each column on the same axes
    for column, color in zip(df.columns, colors):
        hist, bins = np.histogram(df[column], bins=100)
        plt.plot(bins[:-1], hist, color=color, drawstyle='steps-post')
    plt.tight_layout()
    plt.show()

seed = int(input("Please enter the lower value fot the seed: "))
DD = genera_df_3(Quanti_campioni=800, Pronome="25_06_", Size=0.1, Strain_counts=[1,1,1,1,1], Genome_source="Simul", Lower_Seed=seed)

#DD = genera_df_2(1000, "De_caesar_simulated_test_1", ["Aureus_S", "Aureus_R", "Coli", "Pseudomonas", "Arcobacter"], Size = 0.01, Strain_counts = [1,1,1,1,1], Seed=1)
#DD.
DD.to_csv("Da_aggiungere.csv")