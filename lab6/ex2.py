import random, io, textwrap
from Bio import Entrez, SeqIO
import numpy as np
import matplotlib.pyplot as plt

EMAIL = "bugamihai157@gmail.com"
N_GENOMES = 10
RANDOM_SEED = 11

Entrez.email = EMAIL
random.seed(RANDOM_SEED)

def esearch_ids(queries, need):
    seen = []
    for q in queries:
        h = Entrez.esearch(db="nucleotide", term=q, retmax=500)
        ids = Entrez.read(h)["IdList"]
        for i in ids:
            if i not in seen:
                seen.append(i)
            if len(seen) >= need:
                return seen[:need]
    return seen

def fetch_fasta_concat(nid):
    h = Entrez.efetch(db="nucleotide", id=nid, rettype="fasta", retmode="text")
    recs = list(SeqIO.parse(io.StringIO(h.read()), "fasta"))
    if not recs:
        return nid, ""
    seq = "".join(str(r.seq).upper() for r in recs) if len(recs) > 1 else str(recs[0].seq).upper()
    label = recs[0].id
    return label, seq

def ecoRI_fragments(seq):
    motif = "GAATTC"
    cut_offset = 1
    s, sites = 0, []
    while True:
        i = seq.find(motif, s)
        if i == -1:
            break
        sites.append(i + cut_offset)
        s = i + 1
    cuts = [0] + sorted(sites) + [len(seq)]
    return [cuts[i+1] - cuts[i] for i in range(len(cuts)-1)]

def y_from_bp(bp, a, b, lo, hi, H):
    y = a - b * (np.log10(bp) - lo) / (hi - lo + 1e-9)
    return float(np.clip(y, 40, H - 40))

LANE_WIDTH = 20
FIG_DPI = 200
LANE_SPACING = 180
BG_COLOR = "white"
PANEL_COLOR = "#f7f7f7"
LABEL_COLOR = "black"
TICK_COLOR = "black"
TICK_FONT = 18
LANE_LABEL_FONT = 14
SINGLE_TITLE_FONT = 22
ANNOT_LABEL_FONT = 10

def plot_ladder(ax, x_left, a, b, lo, hi, H, fs_tick=TICK_FONT):
    ticks = [12000, 10000, 8000, 6000, 4000, 3000, 2000, 1500, 1000, 700, 500, 300, 150]
    for bp in ticks:
        y = y_from_bp(bp, a, b, lo, hi, H)
        ax.plot([x_left, x_left + 28], [y, y], color=TICK_COLOR, lw=2.0)
        ax.text(x_left - 10, y, f"{bp} bp", va="center", ha="right", color=LABEL_COLOR, fontsize=fs_tick, weight="bold")

def plot_lane(ax, x_center, H, sizes, a, b, lo, hi, color, annotate=False, fsize=ANNOT_LABEL_FONT):
    ax.add_patch(plt.Rectangle((x_center - LANE_WIDTH/2 - 2, 16), LANE_WIDTH + 4, H-32, color=PANEL_COLOR, zorder=0))
    label_y_prev, min_gap = -1e9, 28.0
    for s in sorted(sizes):
        y = y_from_bp(s, a, b, lo, hi, H)
        bh = 8 + max(0, 1000 - s) / 350.0
        ax.add_patch(plt.Rectangle((x_center - LANE_WIDTH/2, y - bh/2), LANE_WIDTH, bh, color=color, zorder=1))
        if annotate:
            y_lab = y if y - label_y_prev >= min_gap else label_y_prev + min_gap
            ax.plot([x_center + LANE_WIDTH/2, x_center + LANE_WIDTH/2 + 12], [y, y_lab], color=LABEL_COLOR, lw=1.5)
            ax.text(x_center + LANE_WIDTH/2 + 14, y_lab, f"{int(s)}", va="center", ha="left",
                    color=LABEL_COLOR, fontsize=fsize,
                    bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="#cccccc", linewidth=0.6))
            label_y_prev = y_lab

queries = [
    'Influenza A virus[ORGN] AND 13000:16000[SLEN]',
    '(influenza[ALL]) AND 12000:17000[SLEN]',
    'Influenza B virus[ORGN] AND 14000:16000[SLEN]',
    'Influenza C virus[ORGN] AND 12000:16000[SLEN]'
]
ids = esearch_ids(queries, N_GENOMES)
if len(ids) < N_GENOMES:
    raise RuntimeError(f"Found {len(ids)} genomes; widen query or check network.")

labels, genomes = [], []
for nid in ids:
    lab, seq = fetch_fasta_concat(nid)
    if seq:
        labels.append(lab)
        genomes.append(seq)
if len(labels) < N_GENOMES:
    raise RuntimeError("Fewer than 10 genomes after fetch; rerun or widen query.")

digests = [ecoRI_fragments(seq) for seq in genomes]
gen_lengths = [len(g) for g in genomes]
all_frags = [f for sizes in digests for f in sizes]

H = 1000
W = 240 + (N_GENOMES+1)*LANE_SPACING
a, b = H * 0.865, H * 0.64
lo = np.log10(max(120, min(all_frags + [150])))
hi = np.log10(max(all_frags + [12000]))

plt.figure(figsize=(W/80, H/80), dpi=FIG_DPI)
ax = plt.gca()
ax.add_patch(plt.Rectangle((0, 0), W, H, color=BG_COLOR))
plot_ladder(ax, 90, a, b, lo, hi, H, fs_tick=TICK_FONT)
for i, (lab, sizes) in enumerate(zip(labels, digests)):
    xc = 200 + i*LANE_SPACING
    color = plt.cm.tab10(i % 10)
    plot_lane(ax, xc, H, sizes, a, b, lo, hi, color=color, annotate=False, fsize=10)
    wrapped = "\n".join(textwrap.wrap(lab, width=20))
    ax.text(xc, H - 28, wrapped, color=LABEL_COLOR, ha="center", fontsize=LANE_LABEL_FONT, rotation=20, weight="semibold")
ax.set_xlim(0, W)
ax.set_ylim(H, 0)
ax.axis("off")
plt.tight_layout()
plt.savefig("gel_combined.png", dpi=FIG_DPI, bbox_inches="tight", facecolor=BG_COLOR)
plt.show()

for i, (lab, sizes, total_len) in enumerate(zip(labels, digests, gen_lengths), start=1):
    h_i, w_i = 1100, 900
    a_i, b_i = h_i*0.865, h_i*0.64
    lo_i = np.log10(max(120, min(sizes + [150])))
    hi_i = np.log10(max(sizes + [12000]))
    plt.figure(figsize=(w_i/80, h_i/80), dpi=FIG_DPI)
    ax = plt.gca()
    ax.add_patch(plt.Rectangle((0, 0), w_i, h_i, color=BG_COLOR))
    plot_ladder(ax, 100, a_i, b_i, lo_i, hi_i, h_i, fs_tick=22)
    plot_lane(ax, 320, h_i, sizes, a_i, b_i, lo_i, hi_i, color="#2b6ea3", annotate=True, fsize=ANNOT_LABEL_FONT)
    title_label = f"{lab}  |  total: {total_len} nt"
    wrapped_title = "\n".join(textwrap.wrap(title_label, width=40))
    ax.text(w_i / 2, h_i - 26, wrapped_title, color=LABEL_COLOR, ha="center", fontsize=SINGLE_TITLE_FONT, weight="bold")
    ax.set_xlim(0, w_i)
    ax.set_ylim(h_i, 0)
    ax.axis("off")
    plt.tight_layout()
    plt.savefig(f"gel_{i}_{lab}.png", dpi=FIG_DPI, bbox_inches="tight", facecolor=BG_COLOR)
    plt.show()

longest_idx = int(np.argmax(gen_lengths))
print("Genome_labels:", labels)
print("Genome_lengths_nt:", gen_lengths)
print("Longest_genome_index_1based:", longest_idx+1)
print("Longest_genome_label:", labels[longest_idx])
print("Longest_genome_length_nt:", gen_lengths[longest_idx])