import numpy as np

def global_alignment(v, w, scoring_matrix, indel_penalty, band_width):
    s = np.zeros((len(v)+1, len(w)+1), dtype=np.int)
    backtrack = {}

    # Initialise first row and column
    for i in range(1, len(v)+1):
        s[i, 0] = s[i-1, 0] - indel_penalty
        backtrack[(i, 0)] = (i-1, 0)
    
    for j in range(1, len(w)+1):
        s[0, j] = s[0, j-1] - indel_penalty
        backtrack[(0, j)] = (0, j-1)

    # Dynamic programming
    for i in range(1, len(v)+1):
        for j in range(max(1, i-band_width), min(len(w)+1, i+band_width+1)):
            # Going down or right
            s[i, j] = max(s[i-1, j] - indel_penalty, s[i, j-1] - indel_penalty)

            # If match or mismatch
            score = scoring_matrix[(v[i-1], w[j-1])]
            s[i, j] = max(s[i, j], \
                          s[i-1, j-1] + score)
    
            # Backtrack to previous coordinate
            if s[i, j] == s[i-1, j-1] + score:
                backtrack[(i, j)] = (i-1, j-1)
            elif s[i, j] == s[i-1, j] - indel_penalty:
                backtrack[(i, j)] = (i-1, j)
            elif s[i, j] == s[i, j-1] - indel_penalty:
                backtrack[(i, j)] = (i, j-1)

    print(s)            
    score = s[len(v), len(w)]

    # Backtrack to find alignments
    i, j = len(v), len(w)
    v_aligned, w_aligned = '', ''
    while i > 0 or j > 0:
        if i > 0 and backtrack[(i, j)] == (i-1, j):
            i-=1
            v_aligned = v[i] + v_aligned
            w_aligned = '-' + w_aligned
        elif j > 0 and backtrack[(i, j)] == (i, j-1):
            j-=1
            v_aligned = '-' + v_aligned
            w_aligned = w[j] + w_aligned
        else:
            i, j = i-1, j-1
            v_aligned = v[i] + v_aligned
            w_aligned = w[j] + w_aligned
    
    return (score, v_aligned, w_aligned)

# Parse the input
f = open('ba5e_small.txt', 'r')
v = f.readline().strip()
w = f.readline().strip()
f.close()

# Parse the scoring matrix
# f = open('blosum62.txt', 'r')
f = open('dna_alignment.txt', 'r')
aas = f.readline().split() # amino acid codes
scoring_matrix = {}

for i in range(len(aas)):
    scores_i = f.readline().split()
    for j in range(len(scores_i)):
        scoring_matrix[(aas[i], aas[j])] = int(scores_i[j])
f.close()


indel_penalty = 1
band_width = 3
(score, v_aligned, w_aligned) = global_alignment(v, w, scoring_matrix, indel_penalty, band_width)
print(score)
print(v_aligned)
print(w_aligned)