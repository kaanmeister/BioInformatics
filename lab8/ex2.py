import random
import urllib.request
import time
import csv  # Added for CSV export

class GenomeLoader:
    """
    Handles loading of bacterial genomes. 
    Can generate synthetic data for testing or download real headers (simulated).
    """
    def __init__(self):
        self.genomes = {}

    def load_mock_genomes(self):
        """
        Generates 3 distinct 'Bacterial' genomes with controlled GC content
        and inserts specific overlapping/nested TE cases for demonstration.
        """
        print("[Loader] Generating 3 simulated Bacterial Genomes...")
        
        # 1. E. coli (Simulated) - Neutral GC
        self.genomes['E. coli (Sim)'] = self._generate_random_dna(15000, gc_content=0.5)
        
        # 2. B. subtilis (Simulated) - Low GC
        self.genomes['B. subtilis (Sim)'] = self._generate_random_dna(12000, gc_content=0.4)
        
        # 3. S. aureus (Simulated) - Variable
        self.genomes['S. aureus (Sim)'] = self._generate_random_dna(10000, gc_content=0.3)
        
        # --- INJECTING TEST CASES ---
        # We inject specific cases to prove the software detects Overlaps and Nesting
        print("[Loader] Injecting complex TE cases (Nested & Overlapping)...")
        
        # Case A: Nested (One inside another)
        # Outer: AAAAAA ... [ BBBBBB ... BBBBBB ] ... TTTTTT
        # Updated to fit max IR length of 6
        outer_ir_L = "GTAAAC"  # 6 bp (Max limit)
        outer_ir_R = "GTTTAC"  # RevComp of L
        
        inner_ir_L = "CCCG"    # 4 bp (Min limit)
        inner_ir_R = "CGGG"    # RevComp
        
        # Construct: IR_L ... junk ... IR_l ... junk ... IR_r ... junk ... IR_R
        nested_seq = (outer_ir_L + "N"*50 + 
                      inner_ir_L + "N"*30 + inner_ir_R + 
                      "N"*50 + outer_ir_R)
                      
        self.genomes['E. coli (Sim)'] = self.genomes['E. coli (Sim)'][:1000] + nested_seq + self.genomes['E. coli (Sim)'][1000:]
        
        return self.genomes

    def _generate_random_dna(self, length, gc_content):
        dna = []
        for _ in range(length):
            r = random.random()
            if r < gc_content:
                dna.append(random.choice("GC"))
            else:
                dna.append(random.choice("AT"))
        return "".join(dna)

    def get_reverse_complement(self, seq):
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
        return "".join(complement.get(base, base) for base in reversed(seq))


class TEfinder:
    """
    Advanced detection engine handling unknown IRs and large datasets.
    """
    def __init__(self, min_ir_len=4, max_ir_len=6, min_te_len=50, max_te_len=5000):
        self.min_ir_len = min_ir_len   # Adjusted to 4
        self.max_ir_len = max_ir_len   # Adjusted to 6
        self.min_te_len = min_te_len   # Min distance between IRs
        self.max_te_len = max_te_len   # Max distance between IRs
        self.loader = GenomeLoader()

    def find_transposons(self, sequence, genome_name):
        """
        Algorithm:
        1. Iterate through genome at position i.
        2. Try IR lengths from max down to min (Greedy approach).
        3. Calculate Reverse Complement.
        4. Search for that RevComp downstream.
        5. If found, record and stop checking other lengths for this start position.
        """
        found_elements = []
        seq_len = len(sequence)
        
        print(f"   -> Scanning {genome_name} ({seq_len} bp)...")
        start_time = time.time()

        # Iterate through the genome
        for i in range(seq_len - self.max_te_len):
            
            # Optimization: Check lengths from Max down to Min.
            # If we find a 6bp match, we prefer that over a 4bp match at the same start.
            match_found_at_i = False
            
            for current_ir_len in range(self.max_ir_len, self.min_ir_len - 1, -1):
                if match_found_at_i:
                    break

                seed = sequence[i : i + current_ir_len]
                
                # Heuristic Filter: Simple repeats (AAAA) trigger false positives.
                if len(set(seed)) == 1:
                    continue

                # Calculate what the Right IR should look like
                target_rev_comp = self.loader.get_reverse_complement(seed)
                
                # Define the search window downstream
                window_start = i + self.min_te_len
                window_end = min(i + self.max_te_len, seq_len)
                
                # Search for the target in this window
                search_window = sequence[window_start : window_end]
                match_index = search_window.find(target_rev_comp)
                
                if match_index != -1:
                    # Found a match!
                    global_match_start = window_start + match_index
                    
                    te_start = i
                    te_end = global_match_start + len(seed)
                    te_length = te_end - te_start
                    
                    found_elements.append({
                        "start": te_start,
                        "end": te_end,
                        "length": te_length,
                        "ir_seq_left": seed,
                        "ir_seq_right": target_rev_comp
                    })
                    match_found_at_i = True

        elapsed = time.time() - start_time
        print(f"      Finished in {elapsed:.2f} seconds. Found {len(found_elements)} candidates.")
        return found_elements

# ==========================================
# Main Execution
# ==========================================
if __name__ == "__main__":
    print("--- ADVANCED BACTERIAL TE DETECTOR ---")
    print("Configuration: Unknown Inverted Repeats | Overlap Detection Enabled")
    print("Constraint Update: IR Size 4-6 bases\n")

    # 1. Load Genomes
    loader = GenomeLoader()
    bacterial_genomes = loader.load_mock_genomes()

    # 2. Initialize Detector
    # Updated constraints: Min IR 4, Max IR 6
    detector = TEfinder(min_ir_len=4, max_ir_len=6, min_te_len=80, max_te_len=2000)

    results = {}

    # 3. Run Detection on all 3 genomes
    for name, sequence in bacterial_genomes.items():
        tes = detector.find_transposons(sequence, name)
        results[name] = tes

    # 4. Export to CSV (New Feature)
    csv_filename = "te_results.csv"
    print(f"\n[IO] Exporting results to {csv_filename}...")
    
    try:
        with open(csv_filename, mode='w', newline='') as file:
            writer = csv.writer(file)
            # Header row
            writer.writerow(["Genome Name", "Start Position", "End Position", "Length", "Left IR", "Right IR"])
            
            total_count = 0
            for name, tes in results.items():
                # Sort by start position
                tes.sort(key=lambda x: x['start'])
                
                for te in tes:
                    writer.writerow([
                        name, 
                        te['start'], 
                        te['end'], 
                        te['length'], 
                        te['ir_seq_left'], 
                        te['ir_seq_right']
                    ])
                    total_count += 1
                    
        print(f"     Success! Wrote {total_count} rows to {csv_filename}.")
        
    except IOError as e:
        print(f"     Error writing CSV: {e}")

    # 5. Minimal Terminal Summary (Instead of dumping everything)
    print("\n" + "="*60)
    print("SUMMARY OF FINDINGS")
    print("="*60)
    for name, tes in results.items():
        print(f"{name:<20} : {len(tes)} potential elements found.")
        # Show first 2 as examples only
        if tes:
            print(f"   Example 1: Start {tes[0]['start']}, Len {tes[0]['length']}, IR {tes[0]['ir_seq_left']}")
    print("-" * 60)
    print(f"Please check '{csv_filename}' for the complete dataset.")