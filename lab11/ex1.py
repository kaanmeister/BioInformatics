import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def needleman_wunsch_with_matrix_and_path(seq1, seq2, match=1, mismatch=-1, gap=0):
    n = len(seq1)
    m = len(seq2)
    
    #the matrix dimensions
    score_matrix = np.zeros((n + 1, m + 1))
    
    # Initialize the first row and column
    for i in range(n + 1):
        score_matrix[i][0] = i * gap
    for j in range(m + 1):
        score_matrix[0][j] = j * gap
        
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if seq1[i - 1] == seq2[j - 1]:
                diagonal_score = score_matrix[i - 1][j - 1] + match
            else:
                diagonal_score = score_matrix[i - 1][j - 1] + mismatch
            
            up_score = score_matrix[i - 1][j] + gap
            left_score = score_matrix[i][j - 1] + gap
            
            score_matrix[i][j] = max(diagonal_score, up_score, left_score)

    align1 = ""
    align2 = ""
    matches = 0
    path = [] 
    
    i, j = n, m
    path.append((i, j))
    
    while i > 0 and j > 0:
        score = score_matrix[i][j]
        score_diag = score_matrix[i - 1][j - 1]
        score_up = score_matrix[i - 1][j]
        score_left = score_matrix[i][j - 1]
        
        if seq1[i-1] == seq2[j-1]:
            step_score = match
        else:
            step_score = mismatch

        # Determine direction
        if score == score_diag + step_score:
            align1 += seq1[i - 1]
            align2 += seq2[j - 1]
            if seq1[i - 1] == seq2[j - 1]:
                matches += 1
            i -= 1
            j -= 1
        elif score == score_up + gap:
            align1 += seq1[i - 1]
            align2 += "-"
            i -= 1
        elif score == score_left + gap:
            align1 += "-"
            align2 += seq2[j - 1]
            j -= 1
        
        path.append((i, j))

    while i > 0:
        align1 += seq1[i - 1]
        align2 += "-"
        i -= 1
        path.append((i, j))
    while j > 0:
        align1 += "-"
        align2 += seq2[j - 1]
        j -= 1
        path.append((i, j))
    
    align1 = align1[::-1]
    align2 = align2[::-1]
    
    return align1, align2, matches, score_matrix, path

#main
S1 = "ACCGTGAAGCCAATAC"
S2 = "AGCGTGCAGCCAATAC"
gap_penalty = 0
match_score = 1
mismatch_penalty = -1

aligned_s1, aligned_s2, match_count, matrix, path = needleman_wunsch_with_matrix_and_path(S1, S2, match_score, mismatch_penalty, gap_penalty)

# --- 1. Text Output ---
visual_match = ""
for k in range(len(aligned_s1)):
    if aligned_s1[k] == aligned_s2[k]:
        visual_match += "|"
    else:
        visual_match += " "

alignment_length = len(aligned_s1)
similarity_percentage = (match_count / alignment_length) * 100

print("Show Alignment:")
print(f"{aligned_s1}")
print(f"{visual_match}")
print(f"{aligned_s2}")
print(f"\nMatches = {match_count}")
print(f"Length = {alignment_length}")
print(f"Similarity = {int(similarity_percentage)} %")

cols = ["-"] + list(S2)
rows = ["-"] + list(S1)
df_matrix = pd.DataFrame(matrix, index=rows, columns=cols)
print("\nScore Matrix:")
print(df_matrix)


# Plot A: Heatmap
plt.figure(figsize=(10, 8))
plt.imshow(matrix, cmap='magma', interpolation='nearest')
plt.title("Graphic representation of the alignment matrix")
plt.xlabel("Sequence 2")
plt.ylabel("Sequence 1")
plt.xticks(np.arange(len(S2)+1), ["-"] + list(S2))
plt.yticks(np.arange(len(S1)+1), ["-"] + list(S1))
plt.colorbar(label='Score')
plt.savefig('heatmap.png')

traceback_grid = np.zeros_like(matrix)
for r, c in path:
    traceback_grid[r, c] = 1

plt.figure(figsize=(10, 8))
# Custom colormap: Yellow background, Red path
from matplotlib.colors import ListedColormap
cmap = ListedColormap(['#FFFFE0', '#D2222D']) 

plt.imshow(traceback_grid, cmap=cmap, interpolation='nearest', aspect='equal')

ax = plt.gca()
ax.set_xticks(np.arange(-0.5, len(S2)+1, 1), minor=True)
ax.set_yticks(np.arange(-0.5, len(S1)+1, 1), minor=True)
ax.grid(which='minor', color='black', linestyle='-', linewidth=1)
ax.tick_params(which='minor', bottom=False, left=False)

plt.title("Traceback path deviation from optimal alignment")
plt.xlabel("Sequence 2")
plt.ylabel("Sequence 1")
plt.xticks(np.arange(len(S2)+1), ["-"] + list(S2))
plt.yticks(np.arange(len(S1)+1), ["-"] + list(S1))

plt.savefig('traceback_path.png')