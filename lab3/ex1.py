def tm_simple(sequence):
    """
    uses the formula: Tm = 4(G + C) + 2(A + T)
    """
    sequence = sequence.upper()
    a = sequence.count('A')
    t = sequence.count('T')
    g = sequence.count('G')
    c = sequence.count('C')
    tm = 4 * (g + c) + 2 * (a + t)
    return tm

def tm_advanced(sequence, na_conc=0.050):
    """
    Uses the formula: 
    Tm = 81.5 + 16.6 * log10([Na+]) + .41 * (%GC) - 600 / length
    """
    from math import log10
    sequence = sequence.upper()
    length = len(sequence)
    gc = sequence.count('G') + sequence.count('C')
    gc_percent = (gc / length) * 100
    tm = 81.5 + 16.6 * log10(na_conc) + 0.41 * gc_percent - 600 / length
    return tm

dna = input("Enter DNA sequence: ")

print("Tm (simple formula): {:.2f} °C".format(tm_simple(dna)))
print("Tm (advanced formula): {:.2f} °C".format(tm_advanced(dna)))
