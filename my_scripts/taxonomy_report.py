import os
import pandas as pd

vocabolario = {
    "LT906455": "Arcobacter",
    "113934c8a6cc489a_1": "Aureus_R", #ATCC BAA-2312
    "113934c8a6cc489a_2": "Aureus_R",
    "c2c883777a574734_1": "Aureus_S", #ATCC 25923
    "c2c883777a574734_2": "Aureus_S",
    "1dad9b1ccf124c63_1": "Coli", #ATCC 25922
    "1dad9b1ccf124c63_2": "Coli",
    "1dad9b1ccf124c63_3": "Coli",
    "1dad9b1ccf124c63_4": "Coli",
    "1dad9b1ccf124c63_5": "Coli",
    "5cfa15b5df924255_1": "Pseudomonas" #ATCC 27853
}

tax_fold = "../risultati_tassonomia_sui_campioni_Decesare"
Abbondanze_secondo_Decesare = pd.read_csv("Abbondanze_secondo_Decesare.csv",header=0, index_col=0)

tax_samples = os.listdir(tax_fold)
Abbondanze_secondo_Ettore = {}
for f in tax_samples:
	counts = pd.read_csv(os.path.join(tax_fold, f), sep="\t", header=0, index_col=0)
	species_counts = {}
	for i in counts.index:
		if vocabolario[i] in species_counts:
			species_counts[vocabolario[i]] += counts.loc[i, "mapped_reads"]	
		else:
			species_counts[vocabolario[i]] = counts.loc[i, "mapped_reads"]
	Abbondanze_secondo_Ettore[f] = species_counts
Abbondanze_secondo_Ettore = pd.DataFrame(Abbondanze_secondo_Ettore, ).T.fillna(0)

