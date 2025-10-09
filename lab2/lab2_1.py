from itertools import product

# DNA bases
bases = ['A', 'C', 'G', 'T']

#dinucleotide combinations (2-mers)
dinucleotides = [''.join(p) for p in product(bases, repeat=2)]
print("All possible dinucleotides:")
print(dinucleotides)

#all trinucleotide combinations (3-mers)
trinucleotides = [''.join(p) for p in product(bases, repeat=3)]
print("\nAll possible trinucleotides:")
print(trinucleotides)
