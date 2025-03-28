import pandas as pd
import matplotlib.pyplot as plt

def plotta_contigs_with_coverage(scaffolds_file, title = "", path = None):
    Alignement_average_quality = None
    contigs_list = []
    with open(scaffolds_file, 'r') as file:
        for line in file:
            if line.startswith('>'):
                contigs_list.append(line.lstrip('>').strip())
    contigs_info = {c : {"len":int(c.split(" ")[1].split("=")[1]) , "cov":float(c.split(" ")[2].split("=")[1][:-1])} for c in contigs_list}
    contigs_info = pd.DataFrame(contigs_info).T
    contigs_info["color"] = "blue"
    max=0
    fig, ax = plt.subplots()
    ax.scatter(contigs_info['len'], contigs_info['cov'], c=contigs_info["color"])
    ax.set_xlabel('Length')
    ax.set_ylabel('Coverage')
    title = title + " Coverage vs Length for " + str(len(contigs_info)) + " contigs"
    ax.set_title(title)
    ax.set_title(title)
    ax.set_xscale("log")
    #ax.set_yscale("log")
    if path is not None:
        plt.savefig(path)
    else:
        plt.show()
    #return contigs_info

assembly_paths = ["Genomes_assembled/assembly/DCWGS1_S1_L001/assembly.fasta",
                "Genomes_assembled/assembly/DCWGS2_S7_L001/assembly.fasta",
                "Genomes_assembled/assembly/DCWGS3_S13_L001/assembly.fasta",
                "Genomes_assembled/assembly/DCWGS4_S19_L001/assembly.fasta",
                "Genomes_assembled/assembly/DCWGS5_S25_L001/assembly.fasta",
                "Genomes_assembled/assembly/DCWGSEC_S42_L001/assembly.fasta",
                "Genomes_assembled/assembly/DCWGSAR_S48_L001/assembly.fasta",
                "Genomes_assembled/assembly/DCWGSPA_S36_L001/assembly.fasta",
                "Genomes_assembled/assembly/DCWGSSAR_S24_L001/assembly.fasta",
                "Genomes_assembled/assembly/DCWGSSAS_S30_L001/assembly.fasta"]

def plot_contig_lengths_barplot(fasta_file, title="", path=None):
    contig_lengths = []
    contig_names = []
    with open(fasta_file, 'r') as file:
        contig_length = 0
        contig_name = ""
        for line in file:
            if line.startswith('>'):
                if contig_length > 0:
                    contig_lengths.append(contig_length)
                    contig_names.append(contig_name)
                contig_name = line[1:].strip()  # Rimuove il carattere '>' e gli spazi
                contig_length = 0
            else:
                contig_length += len(line.strip())
        if contig_length > 0:
            contig_lengths.append(contig_length)
            contig_names.append(contig_name)
    print(fasta_file.split("/")[-1])
    print(sum(contig_lengths))
    #plt.figure(figsize=(10, 6))
    #plt.bar(contig_names, contig_lengths, color='blue', edgecolor='black')
    #plt.xlabel('Contig')
    #plt.ylabel('Lunghezza')
    #plt.title(title + ' lunghezza delle '+ str(len(contig_lengths)) + ' contig')
    #plt.xticks(rotation=90)  # Ruota le etichette dell'asse x per una migliore leggibilità
    #if path is not None:
    #    plt.savefig(path, bbox_inches='tight')  # Salva l'immagine con i margini stretti
    #else:
    #    plt.show()

consensus_paths = ["Consensus_genomes_assembled/DCWGS1_S1_L001.fasta",
                   "Consensus_genomes_assembled/DCWGS2_S7_L001.fasta",
                   "Consensus_genomes_assembled/DCWGS3_S13_L001.fasta",
                   "Consensus_genomes_assembled/DCWGS4_S19_L001.fasta",
                   "Consensus_genomes_assembled/DCWGS5_S25_L001.fasta",
                   "Consensus_genomes_assembled/DCWGSEC_S42_L001.fasta",
                   "Consensus_genomes_assembled/DCWGSAR_S48_L001.fasta",
                   "Consensus_genomes_assembled/DCWGSPA_S36_L001.fasta",
                   "Consensus_genomes_assembled/DCWGSSAR_S24_L001.fasta",
                   "Consensus_genomes_assembled/DCWGSSAS_S30_L001.fasta"]
for path in consensus_paths:
    plot_contig_lengths_barplot(path, title=path.split("/")[1])