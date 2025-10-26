This folder is dedicated for Bioinformatics lab5 work. This task simulates the principle of shotgun DNA sequencing, a fundamental approach in genomics. The main idea is to break a long DNA molecule into many small, random pieces, sequence those fragments, and attempt to computationally reconstruct the original full sequence using the overlaps among these pieces.

Initially, I tried to work with the provided datasets. However, dataset that I installed gave me errors such as tons of "Null" on the `random_reads.txt`. Then, I manually defined the DNA sequences by multiplying them:

```
dna_sequence = 'ACGTACGTACGTACGTACGTACGT' * 100
```

You can find output on `random_reads.txt`.

Ugur Kaan 1242EA.
```