from Bio import SeqIO
from collections import Counter
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio.Seq import Seq

#count codons in a FASTA genome
def count_codons(fasta_file):
    sequences = []
    for record in SeqIO.parse(fasta_file, "fasta"):  
        sequences.append(str(record.seq).upper())
    seq = "".join(sequences) 
    codons = [seq[i:i+3] for i in range(0, len(seq) - 2, 3) if "N" not in seq[i:i+3]]
    return Counter(codons)


# count codons in each genome
covid_codons = count_codons("covid.fasta")
flu_codons = count_codons("influenza.fasta")

# Convert to DataFrame
df_covid = pd.DataFrame(covid_codons.most_common(10), columns=["Codon", "Frequency"])
df_flu = pd.DataFrame(flu_codons.most_common(10), columns=["Codon", "Frequency"])

# plotting for covid
plt.figure(figsize=(8, 5))
sns.barplot(data=df_covid, x="Codon", y="Frequency", palette="Blues_d")
plt.title("Top 10 Most Frequent Codons - COVID-19")
plt.show()

#plotting
plt.figure(figsize=(8, 5))
sns.barplot(data=df_flu, x="Codon", y="Frequency", palette="Oranges_d")
plt.title("Top 10 Most Frequent Codons - Influenza")
plt.show()

# Compare top codons present in both
common_codons = set(df_covid["Codon"]).intersection(set(df_flu["Codon"]))
print("\nCodons common in both genomes (Top 10 sets):", common_codons)

def get_top_amino_acids(codon_counter):
    seq = "".join([codon * count for codon, count in codon_counter.items()])
    protein = str(Seq(seq).translate())
    aa_counts = Counter(protein)
    return aa_counts.most_common(3)

top3_covid = get_top_amino_acids(covid_codons)
top3_flu = get_top_amino_acids(flu_codons)

print("\nTop 3 amino acids (COVID-19):", top3_covid)
print("Top 3 amino acids (Influenza):", top3_flu)

# Prompt suggestion for AI
def make_prompt(amino_acids):
    names = [aa[0] for aa in amino_acids]
    return f"Suggest foods that are low in these amino acids: {', '.join(names)}"

print("\nExample AI prompt (COVID-19):", make_prompt(top3_covid))
print("Example AI prompt (Influenza):", make_prompt(top3_flu))