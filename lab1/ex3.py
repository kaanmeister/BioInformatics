import tkinter as tk
from tkinter import filedialog, messagebox
from collections import Counter

def alphabet(seq: str):
    return sorted(set(seq))

def rel_freq_pct(seq: str):
    n = len(seq)
    if n == 0:
        return {}
    return {k: (v * 100.0) / n for k, v in Counter(seq).items()}

def parse_fasta(text: str):
    lines = text.strip().splitlines()
    if not lines:
        raise ValueError("Empty FASTA content") #throw an error if there is wrong smth with the FASTA
    header = lines[0]
    if not header.startswith('>'):
        raise ValueError("Invalid FASTA: header must start with '>'")
    seq = ''.join(line.strip() for line in lines[1:])
    return header[1:].strip(), seq

#load the FASTA in order to work on GUI.
def load_fasta():
    path = filedialog.askopenfilename(
        title="Select FASTA file",
        filetypes=[("FASTA files", "*.fasta *.fa *.fna *.faa *.fsa"), ("All files", "*.*")]
    )
    if not path:
        return
    try:
        with open(path, 'r', encoding='utf-8') as f:
            content = f.read()
        header, seq = parse_fasta(content)
        seq = ''.join(ch for ch in seq if not ch.isspace())
        alph = alphabet(seq)
        freqs_pct = rel_freq_pct(seq)
        out = []
        out.append(f"Header: {header}")
        out.append(f"Length: {len(seq)}")
        out.append(f"Alphabet: {alph}")
        for k in sorted(freqs_pct):
            out.append(f"{k}: {freqs_pct[k]:.2f}%")
        output_box.config(state='normal')
        output_box.delete('1.0', tk.END)
        output_box.insert(tk.END, "\n".join(out))
        output_box.config(state='disabled')
    except Exception as e:
        messagebox.showerror("Error", str(e))

#GUI Design, for now its simple, we ll see for the next labs.

root = tk.Tk()
root.title("FASTA Analyzer: Alphabet & Percent Composition")

frm = tk.Frame(root, padx=10, pady=10)
frm.pack(fill='both', expand=True)

btn = tk.Button(frm, text="Open FASTA...", command=load_fasta)
btn.pack(anchor='w')

output_box = tk.Text(frm, height=20, width=80, wrap='word')
output_box.pack(fill='both', expand=True, pady=(8,0))
output_box.config(state='disabled')

root.mainloop()