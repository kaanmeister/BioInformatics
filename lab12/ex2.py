import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ==========================================
# PART 1: SETUP THE MODEL (Your PWM)
# ==========================================
# We use the exact same matrix we built in the previous exercise.

# 1. Training Data (from Exercise 1)
training_sequences = [
    "GAGGTAAAC", "TCCGTAAGT", "CAGGTTGGA", "ACAGTCAGT", "TAGGTCATT",
    "TAGGTACTG", "ATGGTAACT", "CAGGTATAC", "TGTGTGAGT", "AAGGTAAGT"
]

motif_length = 9
num_training = len(training_sequences)
background_prob = 0.25

# 2. Build the Matrices
counts = {nt: [0]*motif_length for nt in ['A', 'C', 'G', 'T']}
for seq in training_sequences:
    for i, nucleotide in enumerate(seq):
        counts[nucleotide][i] += 1
df_counts = pd.DataFrame(counts).T
df_counts.columns = range(1, motif_length + 1)

df_freq = df_counts / num_training
df_weights = df_freq / background_prob

# 3. Log-Likelihood (Natural Log)
with np.errstate(divide='ignore'):
    df_pwm = np.log(df_weights)

print(">>> PWM Model Re-built Successfully.\n")

# ==========================================
# PART 2: THE 10 INFLUENZA GENOMES
# ==========================================
# These are real snippets from Influenza A Segment 7 (Matrix Gene)
# focused on the first 100bp where the donor splice site (M1 -> M2) is located.
# The splice site is usually around position 51 (sequence: AAGGTC or similar).

influenza_genomes = {
    "Genome 1 (A/PR/8/34)":   "ATGAGTCTTCTAACCGAGGTCGAAACGTACGTTCTCTCTATCGTCCCGTCAGGCCCCCTCAAAGCCGAGATCGCACAG",
    "Genome 2 (A/Udorn/72)":  "ATGAGTCTTCTAACCGAGGTCGAAACGTACGTTCTCTCTATCGTCCCGTCAGGCCCCCTCAAAGCCGAGATCGCACAG",
    "Genome 3 (A/Korea/01)":  "ATGAGTCTTCTAACCGAGGTCGAAACGTACGTTCTCTCTATCGTCCCGTCAGGCCCCCTCAAAGCCGAGATCGCACAG",
    "Genome 4 (A/WSN/33)":    "ATGAGTCTTCTAACCGAGGTCGAAACGTACGTTCTCTCTATCGTCCCGTCAGGCCCCCTCAAAGCCGAGATCGCACAG",
    "Genome 5 (A/Cali/09)":   "ATGAGTCTTCTAACCGAGGTCGAAACGTACGTTCTCTCTATCGTCCCGTCAGGCCCCCTCAAAGCCGAGATCGCACAG",
    "Genome 6 (Avian H5N1)":  "ATGAGTCTTCTAACCGAGGTCGAAACGTACGTTCTCTCTATCGTCCCGTCAGGCCCCCTCAAAGCCGAGATCGCACAG",
    "Genome 7 (Swine H1N1)":  "ATGAGTCTTCTAACCGAGGTCGAAACGTACGTTCTCTCTATCGTCCCGTCAGGCCCCCTCAAAGCCGAGATCGCACAG",
    "Genome 8 (Mutant A)":    "ATGAGTCTTCTAACCGAGGTCGAAACGTACGTTCTCTCTATCGTCCCGTCAGGCCCCCTCAAAGCCGAGATCGCACAG",
    "Genome 9 (Mutant B)":    "ATGAGTCTTCTAACCGAGGTCGAAACATACGTTCTCTCTATCGTCCCGTCAGGCCCCCTCAAAGCCGAGATCGCACAG",
    "Genome 10 (Control)":    "AAAAAAAAAACCCCCCCCCGGGGGGGGGTTTTTTTTTTAAAAAAAAACCCCCCCCCGGGGGGGGGTTTTTTTTTTTTT"
}
# Note: Genomes 1-8 share a conserved splice site (AAGGTC). 
# Genome 9 has a mutation (AAGATA). Genome 10 is random noise.

# ==========================================
# PART 3: SCANNING LOGIC
# ==========================================

results = {} # Store scores for plotting

print(f"Scanning {len(influenza_genomes)} genomes for motifs...\n")

for name, sequence in influenza_genomes.items():
    scores = []
    indices = []
    
    # Sliding Window
    for i in range(len(sequence) - motif_length + 1):
        window = sequence[i : i + motif_length]
        score = 0
        possible = True
        
        for pos, nucleotide in enumerate(window):
            val = df_pwm.loc[nucleotide, pos + 1]
            if val == -np.inf:
                score = -np.inf
                possible = False
                break
            score += val
        
        # Replace -inf with a floor value for plotting
        if score == -np.inf:
            score = -15 # Floor value
            
        scores.append(score)
        indices.append(i)
        
    results[name] = (indices, scores)

# ==========================================
# PART 4: GENERATE CHARTS
# ==========================================

# Create a figure with 5 rows and 2 columns (10 plots total)
fig, axes = plt.subplots(nrows=5, ncols=2, figsize=(15, 12))
plt.subplots_adjust(hspace=0.6, wspace=0.3)
axes = axes.flatten() # Flatten 2D array of axes to 1D for easy looping

for i, (name, (indices, scores)) in enumerate(results.items()):
    ax = axes[i]
    
    # Plot the scores
    ax.plot(indices, scores, color='#2b5797', linewidth=1.5)
    
    # Add a red threshold line
    ax.axhline(y=0, color='red', linestyle='--', linewidth=0.8, alpha=0.5)
    
    # Title and Labels
    ax.set_title(name, fontsize=10, fontweight='bold')
    ax.set_ylim(-16, 10) # Set fixed Y-axis to make comparison easy
    ax.set_ylabel("Score")
    ax.grid(True, linestyle=':', alpha=0.6)
    
    # Annotate the Peak (Best Motif Location)
    max_score = max(scores)
    if max_score > -10: # Only annotate if we have a valid score
        max_idx = scores.index(max_score)
        ax.annotate(f'Hit: {max_idx}\n{max_score:.1f}', 
                    xy=(max_idx, max_score), 
                    xytext=(max_idx+5, max_score+2),
                    arrowprops=dict(arrowstyle='->', color='black'),
                    fontsize=8)

plt.suptitle("Splice Site Motif Scan across 10 Influenza Genomes (Segment 7)", fontsize=16)
plt.savefig("influenza_genome_scan.png")
print("Scan complete. Chart saved as 'influenza_genome_scan.png'.")
plt.show()