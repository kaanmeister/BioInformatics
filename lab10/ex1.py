import matplotlib.pyplot as plt

# --- Helper Function for Math ---
def get_window_stats(sub_seq):
    """
    Calculates raw CG% and IC for a single short string (the window).
    """
    length = len(sub_seq)
    
    # CG Calculation
    c = sub_seq.count('C')
    g = sub_seq.count('G')
    cg_percent = ((c + g) / length) * 100
    
    # IC Calculation
    counts = {
        'A': sub_seq.count('A'),
        'C': sub_seq.count('C'),
        'G': sub_seq.count('G'),
        'T': sub_seq.count('T')
    }
    numerator = sum(n * (n - 1) for n in counts.values())
    denominator = length * (length - 1)
    ic_val = (numerator / denominator) * 100
    
    return cg_percent, ic_val

# --- Requirement 3: Function to process CpG content ---
def process_cpg_content(sequence, window_size=30):
    cg_list = []
    num_windows = len(sequence) - window_size + 1
    
    for i in range(num_windows):
        window = sequence[i : i + window_size]
        cg_val, _ = get_window_stats(window)
        cg_list.append(cg_val)
        
    average_cg = sum(cg_list) / len(cg_list)
    return average_cg, cg_list

def process_ic_content(sequence, window_size=30):
    """
    Slides a window of 30bp across the sequence.
    Returns the Average IC (Center of Weight) and the list of all window values.
    """
    ic_list = []
    num_windows = len(sequence) - window_size + 1
    
    for i in range(num_windows):
        window = sequence[i : i + window_size]
        _, ic_val = get_window_stats(window)
        ic_list.append(ic_val)
        
    average_ic = sum(ic_list) / len(ic_list)
    return average_ic, ic_list

# --- Main Execution ---
if __name__ == "__main__":
    # 1. Use the sequence S
    S = "CGGACTGATCTATCTAAAAAAAAAAAAAAAAAAAAAAAAAAACGTAGCATCTATCGATCTATCTAGCGATCTATCTACTACG"
    WINDOW = 30
    
    print(f"Analyzing Sequence Length: {len(S)}")
    print(f"Window Size: {WINDOW}")
    print("-" * 40)

    # 3. Process CpG (Get Average and List)
    avg_cg, cg_values = process_cpg_content(S, WINDOW)
    
    # 4. Process IC (Get Average and List)
    avg_ic, ic_values = process_ic_content(S, WINDOW)

    # Print "Each Value" as requested
    print(f"{'Window':<10} {'Seq Segment (first 10)':<25} {'CG %':<10} {'IC Value'}")
    print("-" * 60)
    for i in range(len(cg_values)):
        snippet = S[i : i + 10] + "..."
        print(f"{i+1:<10} {snippet:<25} {cg_values[i]:.2f}       {ic_values[i]:.2f}")
    
    print("-" * 60)
    
    print(f"\nFinal Center of Weight Results:")
    print(f"CG Content Return Value : {avg_cg:.2f} (Target: 29.27)")
    print(f"IC Content Return Value : {avg_ic:.2f} (Target: 27.53)")

    # 5 & 7. Plot the pattern and center
    plt.figure(figsize=(10, 6))
    
    # Plot pattern trajectory
    plt.plot(cg_values, ic_values, marker='o', linestyle='-', color='blue', alpha=0.4, markersize=3, label='Pattern Trajectory')
    
    # Plot Center of Weight
    plt.scatter([avg_cg], [avg_ic], color='red', s=150, zorder=10, label='Center of Weight')
    plt.text(avg_cg + 1, avg_ic, f"Center\n({avg_cg:.2f}, {avg_ic:.2f})", color='red', fontweight='bold')
    
    plt.title("Promoter Pattern Analysis")
    plt.xlabel("C+G Content (%)")
    plt.ylabel("Kappa Index of Coincidence")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.show()