import matplotlib.pyplot as plt
from math import log10

def tm_simple(seq):
    seq = seq.upper()
    return 4 * (seq.count("G") + seq.count("C")) + 2 * (seq.count("A") + seq.count("T"))

def tm_advanced(seq, na_conc=0.050):
    seq = seq.upper()
    length = len(seq)
    gc = seq.count("G") + seq.count("C")
    gc_percent = (gc / length) * 100
    return 81.5 + 16.6 * log10(na_conc) + 0.41 * gc_percent - 600 / length

def read_fasta(filename):
    with open(filename) as f:
        return "".join([line.strip() for line in f if not line.startswith(">")])

def sliding_window_tms(seq, window=9):
    tm1, tm2, positions = [], [], []
    for i in range(len(seq) - window + 1):
        win = seq[i:i+window]
        tm1.append(tm_simple(win))
        tm2.append(tm_advanced(win))
        positions.append(i + window // 2)
    return positions, tm1, tm2

filename = "sequence.fasta"  
dna_seq = read_fasta(filename)
positions, tm1_list, tm2_list = sliding_window_tms(dna_seq, window=9)

plt.plot(positions, tm1_list, label="Tm Simple")
plt.plot(positions, tm2_list, label="Tm Advanced")
plt.xlabel("Position in Sequence")
plt.ylabel("Melting Temp (°C)")
plt.title("Sliding Window Tm Along DNA Sequence")
plt.legend()
plt.show()
