from itertools import product

S = "TACGTGCGCGCGAGCTATCTACTGACTTACGACTAGTGTAGCTGCATCATCGATCGA"
bases = ['A', 'C', 'G', 'T']

#generate all combinations
dinucleotides = [''.join(p) for p in product(bases, repeat=2)]
trinucleotides = [''.join(p) for p in product(bases, repeat=3)]

# this is a helper func for overlapping count
def count_substring(seq, sub):
    return sum(1 for i in range(len(seq) - len(sub) + 1) if seq[i:i+len(sub)] == sub)

total_dinuc = len(S) - 1
total_trinuc = len(S) - 2

# calculate and show the percentages (bulletpoint2)
print("Dinucleotide Percentages:")
for d in dinucleotides:
    count = count_substring(S, d)
    percent = (count / total_dinuc) * 100
    print(f"{d}: {percent:.2f}%")

print("\nTrinucleotide Percentages:")
for t in trinucleotides:
    count = count_substring(S, t)
    percent = (count / total_trinuc) * 100
    print(f"{t}: {percent:.2f}%")
