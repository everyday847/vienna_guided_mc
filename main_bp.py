from vienna_interface import parse_structured_score, parse_score
from rRNA_23S import seq_23S, struct_23S
from rRNA_16S import seq_16S, struct_16S
import numpy as np
import json
import os

def metropolis(temp, new, old):
    """
    Metropolis criterion implementation
    """
    threshold = np.random.uniform()
    factor = np.exp(-(new-old)/temp)
    return factor < threshold

def parens_to_pairmap(struct)
    pairmap = {}
    parens = []
    for i, c in enumerate(struct_23S):
        if c == '(':
            parens.append(i)
        elif c == ')':
            pairmap[i] = parens.pop()
            pairmap[pairmap[i]] = i

if __name__ == '__main__':

    if len(sys.argv) < 2:
        print("must specify subunit: 16S or 23s")
        quit()

    if sys.argv[1] not in ['16S', '23S']:
        print(sys.argv[1], "malformed")
        quit()

    struct = struct_16S if sys.argv[1] == '16S' else struct_23S
    seq = seq_16S if sys.argv[1] == '16S' else seq_23S
    cache_filename = "cache_bp_16S.json" if sys.argv[1] == '16S' else "cache_bp_23S.json"

    pairmap = parens_to_pairmap(struct)

    temp = 2
    seq_scores = {}
    if os.path.exists(cache_filename):
        with open(cache_filename) as f:
            seq_scores = json.loads(f.read())

    # initial score
    score = parse_structured_score(seq, struct) - parse_score(seq)
    print("Starting", score)

    # run one 150-move trajectory, using the cache if available.
    for _ in range(150):
        # pick a random position, and a new random nucleotide identity
        position = np.random.choice(range(len(seq_23S)))
        new_nt = np.random.choice([c for c in 'ACGU' if c != seq_23S[position]])
        
        # update sequence
        new_seq = seq[:position] + new_nt + seq[position+1:]
        
        # mutate base paired as needed
        if position in pairmap:
            new_nt2 = None
            if new_nt == 'G':
                new_nt2 = np.random.choice(['U', 'C'])
            elif new_nt == 'C':
                new_nt2 = 'G'
            elif new_nt == 'U':
                new_nt2 = np.random.choice(['A', 'G'])
            elif new_nt == 'A':
                new_nt2 = 'U'
            else:
                print("Grave issue: new_nt was not ACGU:", new_nt)
                quit()
            new_seq = seq[:pairmap[position]] + new_nt2 + seq[pairmap[position]+1:]

        if new_seq in seq_scores:
            new_score = seq_scores[new_seq]
        else:
            new_score = parse_structured_score(new_seq, struct) - parse_score(new_seq)
            seq_scores[new_seq] = new_score
        if new_score < score:
            # accepted
            seq = new_seq
            score = new_score
        elif metropolis(temp, new_score, score):
            # accepted thermally
            seq = new_seq
            score = new_score
        else:
            # rejected
            pass

    print("Final", new_seq, new_score)

    with open(cache_filename, 'w') as f:
        f.write(json.dumps(seq_scores))
