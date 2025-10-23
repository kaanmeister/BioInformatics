# RNA genetic code mapping from the chart
genetic_code = {
    'UUU': 'Phe', 'UUC': 'Phe', 'UUA': 'Leu', 'UUG': 'Leu', 'UCU': 'Ser', 'UCC': 'Ser',
    'UCA': 'Ser', 'UCG': 'Ser', 'UAU': 'Tyr', 'UAC': 'Tyr', 'UAA': 'Stop', 'UAG': 'Stop',
    'UGU': 'Cys', 'UGC': 'Cys', 'UGA': 'Stop', 'UGG': 'Trp', 'CUU': 'Leu', 'CUC': 'Leu',
    'CUA': 'Leu', 'CUG': 'Leu', 'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
    'CAU': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln', 'CGU': 'Arg', 'CGC': 'Arg',
    'CGA': 'Arg', 'CGG': 'Arg', 'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile', 'AUG': 'Met',
    'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr', 'AAU': 'Asn', 'AAC': 'Asn',
    'AAA': 'Lys', 'AAG': 'Lys', 'AGU': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
    'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val', 'GCU': 'Ala', 'GCC': 'Ala',
    'GCA': 'Ala', 'GCG': 'Ala', 'GAU': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
    'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly'
}

def translate_rna_to_protein(rna_seq):
    rna_seq = rna_seq.upper().replace(" ", "") #so it doesnt really matter if we write upper letters or na
    protein_seq = []
    for i in range(0, len(rna_seq) - 2, 3):
        codon = rna_seq[i:i+3]
        amino_acid = genetic_code.get(codon, '???')
        if amino_acid == 'Stop':
            break
        protein_seq.append(amino_acid)
    return '-'.join(protein_seq)

if __name__ == "__main__":
    rna_input = input("Enter the RNA coding sequence (e.g., AUGGCCAUG...):\n")
    result = translate_rna_to_protein(rna_input)
    print("Protein sequence:", result)
