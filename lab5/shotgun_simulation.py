import random

def get_random_reads(sequence, num_reads=20, min_len=10, max_len=20):
    print(f"Sampling {num_reads} reads from sequence of length {len(sequence)} ...")
    reads = []
    seq_len = len(sequence)
    for i in range(num_reads):
        read_len = random.randint(min_len, max_len)
        start = random.randint(0, seq_len - read_len)
        reads.append(sequence[start:start + read_len])
        if i < 5:  # show preview for first 5
            print(f"Sample read {i+1}: {reads[-1]} (start={start}, length={read_len})")
    print(f"Total reads sampled: {len(reads)}")
    return reads

if __name__ == "__main__":
    # I defined DNA sequence manually due to some "Null" errors on the dataset. 
    dna_sequence = 'ACGTACGTACGTACGTACGTACGT' * 100  

    print("\n--- Shotgun Simulation, Manual Sequence ---")
    print(f"DNA sequence generated. Length: {len(dna_sequence)}")
    print(f"First 40 bases: {dna_sequence[:40]}")
    print(f"Last  40 bases: {dna_sequence[-40:]}\n")

    # generate random reads
    reads = get_random_reads(dna_sequence, num_reads=20, min_len=10, max_len=20)

    with open("random_reads.txt", "w") as out:
        for read in reads:
            out.write(read + "\n")
    print("\nReads written to 'random_reads.txt'.")
    
    #debug for the terminal
    print("\nPreview of the first 5 reads:")
    for i, read in enumerate(reads[:5], 1):
        print(f"Read {i}: {read} (length={len(read)})")

    print(f"\nTotal reads saved: {len(reads)}")
    print("Script completed successfully.")
