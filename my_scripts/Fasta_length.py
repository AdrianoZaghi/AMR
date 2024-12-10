from Bio import SeqIO

def fasta_length(fasta_file):
    total_length = 0
    for record in SeqIO.parse(fasta_file, "fasta"):
        total_length += len(record.seq)
    return total_length

fasta_files = ["../genomes/"+x+".fasta" for x in ["Arcobacter","Aureus_R","Aureus_S","Coli","Pseudomonas"]]
for fasta_file in fasta_files:
     print(fasta_file)
     length = fasta_length(fasta_file)
     print("Length of the FASTA file:", length)