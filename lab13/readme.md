# Bioinformatics Assignments

This repository contains implementations for three assignments regarding Discrete Linear Dynamical Systems and Markov Chain Analysis.

## Dependencies

To run these scripts, you need Python 3 installed along with the following libraries:

    pip install numpy matplotlib

---

## Assignment 1: Matrix Prediction System

Predicts the state of a vector over 5 discrete steps using an arbitrary square matrix and visualizes the results.

### Execution
Run the following script:

    python ex1.py

### Outputs
The script creates a folder named `assignment_results` containing the following visualizations:

**1. Matrix Heatmap**
![Matrix Heatmap](assignment_results/matrix_heatmap.png)

**2. Prediction Trajectory**
![Prediction Trajectory](assignment_results/prediction_trajectory.png)

---

## Assignment 2: DNA Transition Analysis

Generates a random 50-letter DNA sequence and computes the transition probabilities between nucleotides (A, C, G, T).

### Execution
Run the following script:

    python ex2.py

### Outputs
* **Console:** Prints the generated sequence and the formatted probability matrix.
* **File:** Generates `dna_transitions.json` containing the transition matrix.

---

## Assignment 3: Text and Word Transitions

Analyzes an English text (approx. 300 characters), maps unique words to ASCII symbols, and computes transition probabilities between them.

### Execution
Run the following script:

    python ex3.py

### Outputs
* **Console:** Prints the input text, token count, and a sample of the symbol legend.
* **File:** Generates `word_transitions.json` containing:
    * `legend`: Mapping of ASCII symbols to English words.
    * `matrix`: The calculated transition probabilities between symbols.