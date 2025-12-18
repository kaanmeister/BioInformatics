import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# ==========================================
# PART 1: SETUP & DATA
# ==========================================

# The 10 known motif sequences (Training Data)
sequences = [
    "GAGGTAAAC", "TCCGTAAGT", "CAGGTTGGA", "ACAGTCAGT", "TAGGTCATT",
    "TAGGTACTG", "ATGGTAACT", "CAGGTATAC", "TGTGTGAGT", "AAGGTAAGT"
]

# Sequence S for analysis
seq_S = "CAGGTTGGAAACGTAATCAGCGATTACGCATGACGTAA"

# Parameters
motif_length = 9
num_sequences = len(sequences)
background_prob = 0.25

print(f"--- PARAMETERS ---")
print(f"Total Sequences: {num_sequences}")
print(f"Motif Length: {motif_length}")
print(f"Background Probability: {background_prob}\n")

# ==========================================
# PART 2: MATRIX CALCULATIONS
# ==========================================

# 1. Count Matrix
counts = {nt: [0]*motif_length for nt in ['A', 'C', 'G', 'T']}
for seq in sequences:
    for i, nucleotide in enumerate(seq):
        counts[nucleotide][i] += 1
df_counts = pd.DataFrame(counts).T
df_counts.columns = range(1, motif_length + 1)

# 2. Relative Frequencies Matrix
df_freq = df_counts / num_sequences

# 3. Weight Matrix (Odds)
df_weights = df_freq / background_prob

# 4. Log-Likelihood Matrix (Natural Log)
# We use a context manager to ignore divide-by-zero warnings for log(0)
with np.errstate(divide='ignore'):
    df_log = np.log(df_weights)

# ==========================================
# PART 3: SLIDING WINDOW ANALYSIS
# ==========================================

window_data = []

# Loop through Sequence S
for i in range(len(seq_S) - motif_length + 1):
    window = seq_S[i : i + motif_length]
    score = 0
    possible = True
    
    # Calculate score
    for pos, nucleotide in enumerate(window):
        val = df_log.loc[nucleotide, pos + 1]
        if val == -np.inf:
            score = -np.inf
            possible = False
            break
        score += val
    
    # Store result for DataFrame
    window_data.append({
        'Index': i,
        'Window_Sequence': window,
        'Score': score,
        'Is_Impossible': not possible
    })

# Create a clean DataFrame for the sliding window results
df_sliding = pd.DataFrame(window_data)

# ==========================================
# PART 4: TERMINAL OUTPUT
# ==========================================

print(">>> 1. COUNT MATRIX")
print(df_counts)
print("-" * 50)

print(">>> 2. RELATIVE FREQUENCIES MATRIX")
print(df_freq)
print("-" * 50)

print(">>> 3. WEIGHT MATRIX")
print(df_weights)
print("-" * 50)

print(">>> 4. LOG-LIKELIHOOD MATRIX (Natural Log)")
print(df_log)
print("-" * 50)

print(">>> 5. SLIDING WINDOW SCORES")
# formatting the score to string to handle -inf prettily in print
print(df_sliding.to_string(formatters={'Score': lambda x: f"{x:.4f}"})) 

# Identify best hit
max_score = df_sliding.loc[df_sliding['Score'] != -np.inf, 'Score'].max()
best_row = df_sliding[df_sliding['Score'] == max_score].iloc[0]
print("\n" + "="*35)
print(f"RESULT: Strongest Signal at Index {best_row['Index']}")
print(f"Sequence: {best_row['Window_Sequence']}")
print(f"Score:    {max_score:.4f}")
print("="*35 + "\n")

# ==========================================
# PART 5: SAVE TO CSV FILES
# ==========================================

# Saving matrices to CSV
df_counts.to_csv("1_count_matrix.csv")
df_freq.to_csv("2_frequency_matrix.csv")
df_weights.to_csv("3_weight_matrix.csv")
df_log.to_csv("4_log_likelihood_matrix.csv")

# Saving sliding window results to CSV
df_sliding.to_csv("5_sliding_window_scores.csv", index=False)

print(f"Files saved: count_matrix.csv, frequency_matrix.csv, weight_matrix.csv, log_likelihood_matrix.csv, sliding_window_scores.csv")

# ==========================================
# PART 6: PLOTTING THE OUTPUT
# ==========================================

# Prepare data for plotting
scores = df_sliding['Score'].values
indices = df_sliding['Index'].values

# Handle -inf for plotting (replace with a floor value)
finite_scores = scores[scores != -np.inf]
if len(finite_scores) > 0:
    min_finite = min(finite_scores)
    floor_value = min_finite - 5
else:
    floor_value = -10

plot_scores = [s if s != -np.inf else floor_value for s in scores]

plt.figure(figsize=(12, 6))
plt.plot(indices, plot_scores, marker='o', linestyle='-', color='#2b5797', linewidth=2, label='Score')

# Add "Impossible" zone
plt.axhline(y=floor_value + 0.5, color='red', linestyle='--', alpha=0.5, label='Impossible Match Threshold')

plt.title('Sliding Window Log-Likelihood Scores (Sequence S)', fontsize=14)
plt.xlabel('Start Index', fontsize=12)
plt.ylabel('Log-Likelihood Score', fontsize=12)
plt.xticks(indices)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend()

# Annotate Max Score
if len(finite_scores) > 0:
    idx = int(best_row['Index'])
    sc = best_row['Score']
    plt.annotate(f'Peak Signal\nIndex: {idx}\nScore: {sc:.2f}', 
                 xy=(idx, sc), xytext=(idx + 2, sc - 2),
                 arrowprops=dict(facecolor='black', shrink=0.05))

plt.tight_layout()
plt.savefig('sliding_window_plot.png')
print("Plot saved: sliding_window_plot.png")
plt.show()