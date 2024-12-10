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

def genera_df(Quanti_campioni, Pronome, Size, Strain_counts, Genome_source , reads_len, Lower_Seed = 1):
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
        pop = genera_proporzioni(quali_batteri, lista=possible_bacteria)
        Nome = Pronome + str(i)
        provenienza = Genome_source
        sample = [Nome] + pop + [Size, provenienza, "Campione_simulato", Lower_Seed + i + 1] + Strain_counts + [reads_len]
        data.append(sample)
    df = pd.DataFrame(data, columns=["Nome"] + possible_bacteria +["Size","Real/Simul","Original sample","Seed","Aureus_S_N_strains","Aureus_R_N_strains","Coli_N_strains","Pseudomonas_N_strains","Arcobacter_N_strains", "reads_length"])
    df.set_index("Nome", inplace = True)
    return df

def genera_df_multispecie(Quanti_campioni, Pronome, Size, Strain_counts_same_for_all, Genome_source , reads_len, Lower_Seed = 1):
    """Generate "Quanti_campioni" samples with name "Pronome" + i, with size "Size" and strain counts "Strain_counts" and seed "Seed" specified in the arguments.
    The presence of a strain is chosen randomly, first selecting how many of them, and then selecting them from our genome pool.
    The number of symulated strains is fixed for each of the five strains, and is specified in the "Strain_counts" list.
    Specifying a number of strains different from 0 for a strain not present in the sample, will not be a problem: if the bacteria abundance is selected to be 0, nothing of it will go in to the sample.
    """
    data = []
    possible_bacteria = [
        "Aureus_R.fasta",
        "Aureus_S.fasta",
        "Coli.fasta",
        "Pseudomonas.fasta",
        "Arcobacter.fasta",
        "Bacteroides_1_1_6.fasta",
        "Bacteroides_4_3_47.fasta",
        "Bacteroides_capillosus.fasta",
        "Bacteroides_coprocola.fasta",
        "Bacteroides_dorei.fasta",
        "Bacteroides_fragilis.fasta",
        "Bacteroides_plebeius.fasta",
        "Bacteroides_thetaiotaomicron.fasta",
        "Bacteroides_uniformis.fasta",
        "Bacteroides_vulgatus.fasta",
        "Bifidobacterium_adolescentis.fasta",
        "Bifidobacterium_longum.fasta",
        "Bleuta_hansenii.fasta",
        "Clostridium_asparagiforme.fasta",
        "Clostridium_leptum.fasta",
        "Clostridium_phytofermentans.fasta",
        "Clostridium_saccharolyticum.fasta",
        "Ethanoligenens_harbinense.fasta",
        "Faecalibacterium_prausnitzii.fasta",
        "Holdemania_filiformis.fasta",
        "Marvinbryantia_formatexigens.fasta",
        "Ruminococcaceae_bacterium.fasta",
        "Ruminococcus_albus.fasta",
        "Ruminococcus_torques.fasta",
        "Subdoligranulum_variabile.fasta"]
    for i in range(0,Quanti_campioni):
        quanti_batteri = random.randint(1,len(possible_bacteria))
        quali_batteri = random.sample(possible_bacteria, quanti_batteri)
        pop = genera_proporzioni(quali_batteri, lista=possible_bacteria)
        Nome = Pronome + str(i)
        provenienza = Genome_source
        sample = [Nome] + pop + [Size, provenienza, "Campione_simulato", Lower_Seed + i + 1] + [Strain_counts_same_for_all] + [reads_len]
        data.append(sample)
    df = pd.DataFrame(data, columns=["Nome"] + possible_bacteria +["Size","Real/Simul","Original sample","Seed","N_strains_for_all", "reads_length"])
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
#DD = genera_df(Quanti_campioni=1000, Pronome="Halloween_2024_", Size=0.1, Strain_counts=[1,1,1,1,1], Genome_source="Simul_perfected", reads_len=255,Lower_Seed=seed)
DD = genera_df_multispecie(Quanti_campioni=1000, Pronome="22_11_2024_", Size=0.1, Strain_counts_same_for_all=1, Genome_source="Simul_perfected", reads_len=255,Lower_Seed=seed)
DD.to_csv("Da_aggiungere.csv")