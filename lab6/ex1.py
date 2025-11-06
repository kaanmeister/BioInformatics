import random, io
from Bio import Entrez, SeqIO
import numpy as np
import matplotlib.pyplot as plt

EMAIL = "bugamihai157@gmail.com"
SEARCH_TERM = "Escherichia coli[Organism] AND 1000:3000[SLEN]"
LANE_LABEL = "Random fragments"
N_FRAG = 10
FRAG_MIN = 100
FRAG_MAX = 3000
GEL_HEIGHT = 500
GEL_WIDTH = 270
AGAROSE_PCT = 1.0
RANDOM_SEED = 42

Entrez.email = EMAIL
search = Entrez.esearch(db="nucleotide", term=SEARCH_TERM, retmax=50)
ids = Entrez.read(search)["IdList"]
if not ids:
    raise RuntimeError("No records found.")
random.seed(RANDOM_SEED)
chosen_id = random.choice(ids)

handle = Entrez.efetch(db="nucleotide", id=chosen_id, rettype="fasta", retmode="text")
fasta_text = handle.read()
seq_record = SeqIO.read(io.StringIO(fasta_text), "fasta")
seq_len = len(seq_record.seq)

def sample_fragment(seq_len, length_bp):
    L = min(int(length_bp), seq_len)
    start = random.randint(0, max(0, seq_len - L))
    end = start + L
    return (start, end, L)

length_pool = list(range(FRAG_MIN, min(FRAG_MAX, seq_len) + 1))
chosen_lengths = random.sample(length_pool, k=N_FRAG)

fragments = [sample_fragment(seq_len, L) for L in chosen_lengths]
sizes_bp = np.array(sorted(int(L) for (_, _, L) in fragments), dtype=float)

a = GEL_HEIGHT * 0.85
b = GEL_HEIGHT * (0.55 + 0.15 * AGAROSE_PCT)
log_min = np.log10(min(sizes_bp.min(), 150.0))
log_max = np.log10(max(sizes_bp.max(), 3000.0))

def y_from_bp(bp):
    return float(np.clip(a - b * (np.log10(bp) - log_min) / (log_max - log_min + 1e-9), 40, GEL_HEIGHT - 30))

D = np.array([y_from_bp(L) for L in sizes_bp], dtype=float)

plt.figure(figsize=(3.3, 7.6), dpi=150)
ax = plt.gca()
ax.add_patch(plt.Rectangle((0, 0), GEL_WIDTH, GEL_HEIGHT, color="#0a0a0a"))
lane_x = GEL_WIDTH * 0.5
lane_width = 18 #28
ax.add_patch(plt.Rectangle((lane_x - lane_width/2, 10), lane_width, GEL_HEIGHT-20, color="#111111"))

pairs = sorted(zip(D, sizes_bp), key=lambda x: x[0])
label_y_prev = -1e9
min_gap = 20.0
for y, L in pairs:
    band_height = 6 + max(0, 1000 - L) / 400.0
    ax.add_patch(plt.Rectangle((lane_x - lane_width/2, y - band_height/2),
                               lane_width, band_height, color="#e8f5ff"))
    y_lab = y if y - label_y_prev >= min_gap else label_y_prev + min_gap
    ax.plot([lane_x + lane_width/2, lane_x + lane_width/2 + 6],
            [y, y_lab], color="white", lw=0.9)
    ax.text(lane_x + lane_width/2 + 8, y_lab, f"{int(L)} bp",
            va="center", ha="left", color="white", fontsize=8)
    label_y_prev = y_lab

ladder = [3000, 1500, 1000, 700, 500, 300, 150]
for bp in ladder:
    y = y_from_bp(bp)
    ax.plot([60, 85], [y, y], color="white", lw=1.2, solid_capstyle="butt", clip_on=False)
    ax.text(55, y, f"{bp} bp", va="center", ha="right", color="white", fontsize=10, clip_on=False)

ax.text(lane_x, GEL_HEIGHT - 12, LANE_LABEL, color="white", ha="center", fontsize=10)
ax.set_xlim(0, GEL_WIDTH)
ax.set_ylim(GEL_HEIGHT + 10, -10)
ax.axis("off")
plt.subplots_adjust(left=0.30, right=0.98, top=0.98, bottom=0.08)
plt.tight_layout()
plt.show()

print(f"NCBI ID: {chosen_id}")
print(f"Sequence length: {len(seq_record.seq)} nt")
print("Distinct fragment sizes (bp):", [int(x) for x in sizes_bp])