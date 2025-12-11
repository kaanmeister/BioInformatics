import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch  # <--- FIXED IMPORT
from Bio import Entrez, SeqIO
from scipy.signal import find_peaks
import os

# --- 1. CONFIGURATION ---
Entrez.email = "kaanmeister.ro@gmail.com" 

def download_genome(accession_id, filename):
    if os.path.exists(filename):
        # print(f"Loading {filename}...") 
        record = SeqIO.read(filename, "fasta")
    else:
        print(f"Downloading {accession_id}...")
        try:
            handle = Entrez.efetch(db="nucleotide", id=accession_id, rettype="fasta", retmode="text")
            record = SeqIO.read(handle, "fasta")
            handle.close()
            with open(filename, "w") as output_handle:
                SeqIO.write(record, output_handle, "fasta")
        except Exception as e:
            print(f"Error downloading {accession_id}: {e}")
            return ""
    return str(record.seq)

def calculate_gc_profile(sequence, window_size=300):
    """Calculates GC content sliding window to create a 'profile' for the genome."""
    profile = []
    # Step of 100 for speed and smoother graph
    for i in range(0, len(sequence) - window_size, 100): 
        chunk = sequence[i:i+window_size]
        gc_count = chunk.count('G') + chunk.count('C')
        profile.append(gc_count / window_size)
    return profile

# --- 2. GET DATA ---
# COVID
covid_seq = download_genome("NC_045512", "covid19.fasta")

# INFLUENZA (Concatenate all segments to make one long "genome" for comparison)
flu_accessions = ["NC_002019", "NC_002023", "NC_002021", "NC_002017", 
                  "NC_002022", "NC_002018", "NC_002016", "NC_002020"]
flu_full_seq = ""
segment_boundaries = [] 
current_pos = 0

print("Downloading and assembling Influenza segments...")
for acc in flu_accessions:
    seq = download_genome(acc, f"{acc}.fasta")
    flu_full_seq += seq
    segment_boundaries.append(current_pos)
    current_pos += len(seq)

# --- 3. GENERATE PROFILES ---
print("Calculating alignment profiles...")
covid_profile = calculate_gc_profile(covid_seq)
flu_profile = calculate_gc_profile(flu_full_seq)

# Normalize lengths for plotting side-by-side
x_covid = np.linspace(0, len(covid_seq), len(covid_profile))
x_flu = np.linspace(0, len(flu_full_seq), len(flu_profile))

# --- 4. PLOTTING ---
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=False)

# Plot Influenza (Top)
ax1.plot(x_flu, flu_profile, color='blue', label='Influenza A (Concatenated Segments)')
ax1.set_ylabel("GC Profile")
ax1.set_title("Genome 1: Influenza A (Segments 1-8)")
ax1.set_ylim(0.2, 0.6)
ax1.grid(True, alpha=0.3)

# Draw lines for segments
for boundary in segment_boundaries:
    ax1.axvline(boundary, color='red', linestyle='--', alpha=0.5)
    if boundary > 0: 
        ax1.text(boundary - 50, 0.55, "|", fontsize=8, color='red')

# Plot COVID (Bottom)
ax2.plot(x_covid, covid_profile, color='green', label='SARS-CoV-2')
ax2.set_ylabel("GC Profile")
ax2.set_xlabel("Position (bp)")
ax2.set_title("Genome 2: SARS-CoV-2")
ax2.set_ylim(0.2, 0.6)
ax2.grid(True, alpha=0.3)

# VISUAL ALIGNMENT CONNECTOR
# Find peaks to simulate 'features' to align
peaks_flu, _ = find_peaks(flu_profile, height=0.45, distance=20)
peaks_cov, _ = find_peaks(covid_profile, height=0.40, distance=50)

print("Drawing alignment connections...")

# Draw a few connection lines between features
# (Limiting to 5 lines for clarity)
num_lines = min(5, len(peaks_flu), len(peaks_cov))

for i in range(num_lines):
    flu_x = x_flu[peaks_flu[i]]
    cov_x = x_covid[peaks_cov[i]]
    
    # Use the imported ConnectionPatch
    con = ConnectionPatch(xyA=(flu_x, 0.25), xyB=(cov_x, 0.55), 
                          coordsA="data", coordsB="data", 
                          axesA=ax1, axesB=ax2, color="gray", alpha=0.5, linestyle="-.")
    ax2.add_artist(con)

plt.tight_layout()
plt.savefig("genome_alignment_profile.png")
print("Done! Plot saved as 'genome_alignment_profile.png'")
plt.show()