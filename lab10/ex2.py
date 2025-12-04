import time
import matplotlib.pyplot as plt
from Bio import Entrez, SeqIO
import numpy as np

dataset = {
    # --- 10 COVID-19 Variants (Full Genomes) ---
    "COVID-19 (Wuhan-Hu-1)": "NC_045512",
    "COVID-19 (Alpha)": "MW320722",
    "COVID-19 (Beta)": "MW598419",
    "COVID-19 (Delta)": "MZ558051",
    "COVID-19 (Omicron BA.1)": "OL672836",
    "COVID-19 (Gamma)": "OK173809",
    "COVID-19 (Mu)": "OX315743",
    "COVID-19 (Lambda)": "OV065430",
    "COVID-19 (USA-WA1)": "MN985325",
    "COVID-19 (Italy)": "MT066156",

    # --- 10 Influenza Variants (Hemagglutinin Gene) ---
    "Flu A (1918 H1N1 Spanish)": "AF117241",
    "Flu A (2009 H1N1 Swine)": "CY041645",
    "Flu A (H3N2 Hong Kong)": "CY163680",
    "Flu A (H5N1 Avian)": "AF144305",
    "Flu A (H7N9)": "KC853767",
    "Flu A (H9N2)": "MK516997",
    "Flu B (Yamagata)": "CY115151",
    "Flu B (Victoria)": "CY115159",
    "Flu A (H2N2 1957)": "CY022013",
    "Flu A (H7N7)": "CY009616"
}

def get_window_stats(sub_seq):
    length = len(sub_seq)
    if length == 0: return 0, 0
    
    # CG Calc
    c = sub_seq.count('C')
    g = sub_seq.count('G')
    cg_percent = ((c + g) / length) * 100
    
    # IC Calc
    counts = {'A': sub_seq.count('A'), 'C': sub_seq.count('C'), 
              'G': sub_seq.count('G'), 'T': sub_seq.count('T')}
    numerator = sum(n * (n - 1) for n in counts.values())
    denominator = length * (length - 1)
    ic_val = (numerator / denominator) * 100 if denominator > 0 else 0
    
    return cg_percent, ic_val

def analyze_genome(sequence, window_size=200):
    # Note: Increased window size to 200 for full genomes to reduce noise in plots
    cg_vals = []
    ic_vals = []
    
    seq_str = str(sequence).upper()
    seq_len = len(seq_str)
    
    # Step through sequence
    step = 50 # Optimization: Step size to speed up processing of large genomes
    for i in range(0, seq_len - window_size, step):
        window = seq_str[i : i + window_size]
        cg, ic = get_window_stats(window)
        cg_vals.append(cg)
        ic_vals.append(ic)
        
    avg_cg = sum(cg_vals) / len(cg_vals) if cg_vals else 0
    avg_ic = sum(ic_vals) / len(ic_vals) if ic_vals else 0
    
    return cg_vals, ic_vals, avg_cg, avg_ic

# --- 3. MAIN EXECUTION ---

def run_assignment():
    print("--- Starting Bioinformatics Analysis ---")
    print("Downloading sequences from NCBI... (This might take a minute)")
    
    results = {}
    
    for name, accession in dataset.items():
        try:
            print(f"Fetching {name} ({accession})...")
            # Fetch from NCBI
            handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
            record = SeqIO.read(handle, "fasta")
            handle.close()
            
            # Analyze
            cg_trace, ic_trace, avg_cg, avg_ic = analyze_genome(record.seq)
            
            results[name] = {
                "cg_trace": cg_trace,
                "ic_trace": ic_trace,
                "center_cg": avg_cg,
                "center_ic": avg_ic,
                "type": "Covid" if "COVID" in name else "Flu"
            }
            
        except Exception as e:
            print(f"Error fetching {name}: {e}")

    # --- PLOTTING ---
    
    plt.figure(figsize=(12, 6))
    for name, data in results.items():
        color = 'red' if data['type'] == 'Covid' else 'blue'
        alpha = 0.1 # Make lines very faint so we can see the 'cloud' shape
        plt.plot(data['cg_trace'], data['ic_trace'], color=color, alpha=alpha, linewidth=1)
    
    from matplotlib.lines import Line2D
    custom_lines = [Line2D([0], [0], color='red', lw=2),
                    Line2D([0], [0], color='blue', lw=2)]
    plt.legend(custom_lines, ['COVID-19 Genomes', 'Influenza Genomes'])
    
    plt.title("Chart 1: Objective Digital Stains (Trajectories)")
    plt.xlabel("C+G Content (%)")
    plt.ylabel("Kappa Index of Coincidence (IC)")
    plt.grid(True, alpha=0.3)
    plt.show()

    plt.figure(figsize=(12, 8))
    
    for name, data in results.items():
        color = 'red' if data['type'] == 'Covid' else 'blue'
        marker = 'o' if data['type'] == 'Covid' else '^'
        
        plt.scatter(data['center_cg'], data['center_ic'], color=color, marker=marker, s=100)
        
        # Add labels
        offset_y = 0.2 if data['type'] == 'Covid' else -0.2
        plt.text(data['center_cg'], data['center_ic'] + offset_y, name, 
                 fontsize=8, ha='center', va='center')

    plt.title("Chart 2: Center of Weight Comparison")
    plt.xlabel("Average C+G Content (%)")
    plt.ylabel("Average Kappa IC")
    plt.grid(True, linestyle='--')
    
    plt.scatter([], [], c='red', marker='o', label='COVID-19 Variants')
    plt.scatter([], [], c='blue', marker='^', label='Influenza Variants')
    plt.legend()
    
    plt.show()

if __name__ == "__main__":
    run_assignment()