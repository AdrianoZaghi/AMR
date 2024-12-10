# %%
import M_learning_images as im
import pandas as pd
import matplotlib.pyplot as plt

Campioni = pd.read_csv("Database_campioni2.csv", index_col=0)
campioni = Campioni[(Campioni["Aureus_S_N_strains"] == 1) | (Campioni["Aureus_S_N_strains"] == 0)]
campioni = campioni[(campioni["Aureus_R_N_strains"] == 1) | (campioni["Aureus_R_N_strains"] == 0)]
campioni = campioni[campioni["Size"] == 0.1]
# %%
Isolate = "Pseudomonas"
t = campioni[campioni[Isolate] >= 0 ]

Simulati_PATRIC = ["samples_2/" + x for x in t[t["Real/Simul"] == "Simul"].index]
Simulati_De_Cesare = ["samples_2/" + x for x in t[t["Real/Simul"] == "De_Caesar_simul"].index]
Simulati_new_error = ["samples_2/" + x for x in t[t["Real/Simul"] == "Simul_new_error"].index]
Simulati_new_trim = ["samples_2/" + x for x in t[t["Real/Simul"] == "Simul_new_trim"].index]
Simulati_perfected = ["samples_2/" + x for x in t[t["Real/Simul"] == "Simul_perfected"].index]


t = campioni[campioni[Isolate] >= 0]
Veri = ["samples_2/" + x for x in t[t["Real/Simul"] == "Real"].index]
Veri_new_trim = ["samples_2/" + x for x in t[t["Real/Simul"] == "Real_new_trim"].index]

pro = ["SAR", "SAS", "SEC", "SPA", "SEC"]
Veri_isolati = [x for x in Veri if any([y in x for y in pro])]
Veri_non_isolati = [x for x in Veri if not any([y in x for y in pro])]

# %%
for i in Simulati_perfected[:4]:
    im.read_bwt_2_aro_MAPQ(i, just_completely_mapped=True, cumulative=True)
# %%
for i in Veri_isolati[:4]:
    im.read_bwt_2_aro_MAPQ(i, just_completely_mapped=True, cumulative=True)
# %%
MPQC = -1
Completely_mapped = True
gene_counts_veri = dict(im.read_bwt_2_aro(Veri, 1, mapq_average_min=MPQC, just_completely_mapped=Completely_mapped).sum(axis=0) / (2*len(Veri)))
gene_counts_simulati_PATRIC = dict(im.read_bwt_2_aro(Simulati_PATRIC, 1, mapq_average_min=MPQC, just_completely_mapped=Completely_mapped).sum(axis=0) / len(Simulati_PATRIC))
gene_counts_simulati_De_Cesare = dict(im.read_bwt_2_aro(Simulati_De_Cesare, 1, mapq_average_min=MPQC, just_completely_mapped=Completely_mapped).sum(axis=0) / len(Simulati_De_Cesare))
gene_counts_new_error = dict(im.read_bwt_2_aro(Simulati_new_error, 1, mapq_average_min=MPQC, just_completely_mapped=Completely_mapped).sum(axis=0) / len(Simulati_new_error))
gene_counts_new_trim = dict(im.read_bwt_2_aro(Simulati_new_trim, 1, mapq_average_min=MPQC, just_completely_mapped=Completely_mapped).sum(axis=0) / len(Simulati_new_trim))
gene_counts_perfected = dict(im.read_bwt_2_aro(Simulati_perfected, 1, mapq_average_min=MPQC, just_completely_mapped=Completely_mapped).sum(axis=0) / len(Simulati_perfected))
gene_counts_veri_new_trim = dict(im.read_bwt_2_aro(Veri_new_trim, 1, mapq_average_min=MPQC, just_completely_mapped=Completely_mapped).sum(axis=0) / (len(Veri_new_trim)*2))

df = pd.DataFrame({
    #'Veri': gene_counts_veri,
    'Veri_new_trim': gene_counts_veri_new_trim,
    #'Simulati_PATRIC': gene_counts_simulati_PATRIC,
    #'Simulati_De_Cesare': gene_counts_simulati_De_Cesare,
    #'Simulati_new_error': gene_counts_new_error,
    'Simulati_perfected': gene_counts_perfected
})
# %%
df = df.fillna(0)
df['row_sum'] = df.sum(axis=1)
print(df.sum(axis=0))
#df = df.div(df.sum(axis=0), axis=1)
# Ordina il DataFrame in base alla somma delle righe
df = df.sort_values(by='Simulati_perfected', ascending=False)
df = df.drop(columns=["row_sum"])

df_head = df.head(40)
ax = df_head.plot(kind='bar', figsize=(10, 6), width=0.8)
# Aggiungi etichette e titolo
plt.title('Geni per ciascun Isolato')
plt.xlabel('Geni')
plt.ylabel('Somma')
plt.xticks(rotation=45)
# Mostra la legenda
plt.legend(title='Gruppi')
# Mostra il grafico
plt.show()

# %%
