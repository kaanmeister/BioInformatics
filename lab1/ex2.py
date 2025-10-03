from collections import Counter

def alphabet(seq: str):
    n = len(seq)
    return {k: v / n for k, v in Counter(seq).items()}
    #return sorted(set(seq))

S = "ATTGCCCCGAAT"
print(alphabet(S)) 