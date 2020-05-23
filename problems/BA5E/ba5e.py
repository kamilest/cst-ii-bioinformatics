import numpy as np

def dna_scoring(a, b, match_score, mismatch_penalty):
    if (a == 'A' and b == 'T') \
            or (a == 'T' and b == 'A') \
            or (a == 'C' and b == 'G') \
            or (a == 'G' and b == 'C'):
        return match_score
    else:
        return mismatch_penalty


def global_alignment(v, w, match_score, mismatch_penalty, indel_penalty, scoring_matrix=None):
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
        for j in range(1, len(w)+1):
            # Going down or right
            s[i, j] = max(s[i-1, j] - indel_penalty, s[i, j-1] - indel_penalty)

            # If match or mismatch
            if scoring_matrix is not None:
                score = scoring_matrix[(v[i-1], w[j-1])]
            else: 
                score = dna_scoring(v[i-1], w[j-1], match_score, -1 * mismatch_penalty)
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


def parse_scoring_matrix():
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
    return scoring_matrix

# Parse the input
f = open('ba5e.txt', 'r')
v = f.readline().strip()
w = f.readline().strip()
f.close()

# ALL MUST BE POSITIVE
match_score = 1
mismatch_penalty = 1
indel_penalty = 4

(score, v_aligned, w_aligned) = global_alignment(v, w, match_score, mismatch_penalty, indel_penalty)

print(score)
print(v_aligned)
print(w_aligned)