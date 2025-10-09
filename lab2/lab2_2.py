S = "TACGTGCGCGCGAGCTATCTACTGACTTACGACTAGTGTAGCTGCATCATCGATCGA"

def observed_kmers(seq, k):
    return {seq[i:i+k] for i in range(len(seq)-k+1)}

dinucs = observed_kmers(S, 2)
trinucs = observed_kmers(S, 3)

print("Observed dinucleotides (k=2):")
print(", ".join(sorted(dinucs)))

print("\nObserved trinucleotides (k=3):")
print(", ".join(sorted(trinucs)))
