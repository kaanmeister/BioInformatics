import random

class DNABuilder:
    """
    Handles the generation of random DNA and the biological simulation 
    of Transposable Element (TE) insertion.
    """
    def __init__(self, length=300):
        self.length = length
        self.sequence = ""
        self.ground_truth = [] 
        self.ir_left = "CAGTGC" 
        self.ir_right = "GCACTG" 
        self.te_min_len = 30    
        self.te_max_len = 50    
        self.tsd_len = 4         

    def generate_random_sequence(self, length):
        """Generates a random string of nucleotides."""
        return "".join(random.choice("ATCG") for _ in range(length))

    def build_host(self):
        """Creates the initial host DNA sequence."""
        self.sequence = self.generate_random_sequence(self.length)
        print(f"\n[1] Generated Host DNA ({len(self.sequence)} bp)")
        print(f"    Start: {self.sequence[:10]}... End: ...{self.sequence[-10:]}")

    def insert_transposons(self, num_elements=3):
        """
        Simulates insertion.
        Mechanism: 
        1. pick a random insertion site.
        2. create the TE (Left IR + Random Internal + Right IR).
        3. simulate Target Site Duplication (TSD) (Direct Repeats from your image).
        """
        print(f"\n[2] Inserting {num_elements} Transposable Elements...")
        
        for i in range(num_elements):
            # pick a random insertion point in the CURRENT sequence
            # We avoid the very edges to allow for TSDs
            insert_pos = random.randint(10, len(self.sequence) - 10)
            
            # 2. Create the TE payload
            internal_len = random.randint(self.te_min_len, self.te_max_len)
            internal_seq = self.generate_random_sequence(internal_len)
            te_sequence = self.ir_left + internal_seq + self.ir_right
            
            upstream = self.sequence[:insert_pos]
            downstream = self.sequence[insert_pos:]
            
            # The 'target' that gets duplicated is the few bases AT the insertion point
            target_site = self.sequence[insert_pos : insert_pos + self.tsd_len]
            
            self.sequence = upstream + te_sequence + target_site + downstream[self.tsd_len:]
            
            # Calculate actual start/end for Ground Truth
            # Start is after the first TSD (or including it? usually TE is defined by IRs)
            # Let's define position as the exact start of the TE (IR_Left)
            te_start = insert_pos
            te_end = insert_pos + len(te_sequence)
            
            self.ground_truth.append({
                "id": i+1,
                "start": te_start,
                "end": te_end,
                "seq": te_sequence
            })
            
            print(f"    -> Inserted TE #{i+1} at index {te_start}. Length: {len(te_sequence)}bp")

        self.ground_truth.sort(key=lambda x: x['start'])
        return self.sequence

class TransposonDetector:
    """
    Software to detect the positions of the transposable elements.
    """
    def __init__(self, sequence, ir_left, ir_right):
        self.sequence = sequence
        self.ir_left = ir_left
        self.ir_right = ir_right
    
    def run(self):
        """
        Scans sequence for the Inverted Repeat signature.
        """
        detected_elements = []
        seq_len = len(self.sequence)
        ir_len = len(self.ir_left)
        
        i = 0
        while i < seq_len - ir_len:
            # Check if we hit a Left IR
            chunk = self.sequence[i : i + ir_len]
            
            if chunk == self.ir_left:
                search_limit = min(i + 150, seq_len)
                found_pair = False
                
                # Scan ahead for the Right IR
                for j in range(i + ir_len, search_limit):
                    chunk_right = self.sequence[j : j + ir_len]
                    if chunk_right == self.ir_right:
                        # FOUND ONE!
                        te_start = i
                        te_end = j + ir_len
                        detected_elements.append((te_start, te_end))
                        
                        #move our main pointer past this element to avoid overlapping detection
                        i = te_end
                        found_pair = True
                        break
                
                if not found_pair:
                    i += 1
            else:
                i += 1
                
        return detected_elements

# =============
# Main Execution
# ==============
if __name__ == "__main__":
    
    builder = DNABuilder(length=400) # 200-400b as requested
    builder.build_host()
    final_sequence = builder.insert_transposons(num_elements=3)
    
    print(f"\nFinal Sequence Length: {len(final_sequence)} bp")
    print(f"Final Sequence Preview: {final_sequence[:50]}...{final_sequence[-50:]}")
    
    #run detection
    detector = TransposonDetector(final_sequence, builder.ir_left, builder.ir_right)
    results = detector.run()
    
    #compare the results
    print("\n[4] Results Verification")
    print(f"{'ID':<5} {'Type':<10} {'Start':<10} {'End':<10} {'Status'}")
    print("-" * 50)
    
    #print what we detected
    for idx, (start, end) in enumerate(results):
        # Check against ground truth
        match = "UNKNOWN"
        for gt in builder.ground_truth:
            if gt['start'] == start and gt['end'] == end:
                match = "CORRECT"
                break
        print(f"{idx+1:<5} {'Detected':<10} {start:<10} {end:<10} {match}")
        
    #print what was actually there (Ground Truth)
    print("\n--- Ground Truth Data ---")
    for gt in builder.ground_truth:
        print(f"TE #{gt['id']}: Start {gt['start']}, End {gt['end']}, Seq: {gt['seq'][:10]}...")