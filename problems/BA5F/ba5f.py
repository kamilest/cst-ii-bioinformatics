import numpy as np

def dna_scoring(a, b, match_score, mismatch_penalty):
    if (a == 'A' and b == 'T') \
            or (a == 'T' and b == 'A') \
            or (a == 'C' and b == 'G') \
            or (a == 'G' and b == 'C'):
        return match_score
    else:
        return mismatch_penalty

def local_alignment(v, w, match_score, mismatch_penalty, indel_penalty, scoring_matrix=None):
    s = np.zeros((len(v)+1, len(w)+1), dtype=np.int)
    backtrack = {}

    # Initialise backtrack (scores are already 0)
    for i in range(len(v)+1):
        for j in range(len(w)+1):
            backtrack[(i, j)] = (0, 0)
    

    # Dynamic programming    
    for i in range(1, len(v)+1):
        for j in range(1, len(w)+1):

            if scoring_matrix is not None:
                score = scoring_matrix[(v[i-1], w[j-1])]
            else: 
                score = dna_scoring(v[i-1], w[j-1], match_score, -1 * mismatch_penalty)

            s[i, j] = np.amax([0, \
                        s[i-1, j] - indel_penalty, \
                        s[i, j-1] - indel_penalty, \
                        s[i-1, j-1] + score])
    
            # Backtrack to previous coordinate
            if s[i, j] == s[i-1, j-1] + score:
                backtrack[(i, j)] = (i-1, j-1)
            elif s[i, j] == s[i-1, j] - indel_penalty:
                backtrack[(i, j)] = (i-1, j)
            elif s[i, j] == s[i, j-1] - indel_penalty:
                backtrack[(i, j)] = (i, j-1)

    print(s)

    # Find optimal score and one of the optimal coords.
    score_opt = np.amax(s)
    coords_opt = np.where(s == score_opt)
    (i_opt, j_opt) = list(zip(coords_opt[0], coords_opt[1]))[0]

    # Backtrack to find alignments
    i, j = i_opt, j_opt
    v_aligned, w_aligned = '', ''
    while i > 0 and j > 0:
        if backtrack[(i, j)] == (0, 0):
            i, j = 0, 0
        elif backtrack[(i, j)] == (i-1, j):
            i-=1
            v_aligned = v[i] + v_aligned
            w_aligned = '-' + w_aligned
        elif backtrack[(i, j)] == (i, j-1):
            j-=1
            v_aligned = '-' + v_aligned
            w_aligned = w[j] + w_aligned
        elif backtrack[(i, j)] == (i-1, j-1):
            i, j = i-1, j-1
            v_aligned = v[i] + v_aligned
            w_aligned = w[j] + w_aligned
    
    return (int(score_opt), v_aligned, w_aligned)


def parse_input_strings():
    # Parse the input
    f = open('ba5e.txt', 'r')
    v = f.readline().strip()
    w = f.readline().strip()
    f.close()
    return v, w


def parse_scoring_matrix():
    # Parse the scoring matrix
    f = open('pam250.txt', 'r')
    aas = f.readline().split() # amino acid codes
    scoring_matrix = {}

    for i in range(len(aas)):
        scores_i = f.readline().split()
        for j in range(len(scores_i)):
            scoring_matrix[(aas[i], aas[j])] = int(scores_i[j])
        f.close()

    return scoring_matrix

# ALL MUST BE POSITIVE
match_score = 1
mismatch_penalty = 1
indel_penalty = 4

v = 'GCACTT'
w = 'CCCAAT'

(score, v_aligned, w_aligned) = local_alignment(v, w, match_score, mismatch_penalty, indel_penalty)
print(score)
print(v_aligned)
print(w_aligned)