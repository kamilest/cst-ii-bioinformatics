import numpy as np

def overlap_alignment(v, w, match_score=1, mismatch_penalty=2, indel_penalty=2):
    s = np.zeros((len(v)+1, len(w)+1))
    backtrack = {}

    # Initialise backtrack, scores are already 0
    for i in range(1, len(v)+1):
        for j in range(1, len(w)+1):
            backtrack[(i, j)] = (0, 0)

    # Dynamic programming    
    for i in range(1, len(v)+1):
        for j in range(1, len(w)+1):
            # Match or mismatch
            score = match_score if v[i-1] == w[j-1] else -1 * mismatch_penalty

            # Indels
            s[i, j] = np.amax([s[i-1, j] - indel_penalty, \
                               s[i, j-1] - indel_penalty, \
                               s[i-1, j-1] + score])
            
    
            # Backtrack to previous coordinate
            if s[i, j] == s[i-1, j-1] + score:
                backtrack[(i, j)] = (i-1, j-1)
            elif s[i, j] == s[i-1, j] - 2:
                backtrack[(i, j)] = (i-1, j)
            elif s[i, j] == s[i, j-1] - 2:
                backtrack[(i, j)] = (i, j-1)

                
    # Find optimal score and one of the optimal coords.
    score_opt = max(np.amax(s[:, len(w)]),np.amax(s[len(v), :]))

    if score_opt == np.amax(s[:, len(w)]):
        coords_opt = np.where(s[:, len(w)] == score_opt)
        i_opt, j_opt = coords_opt[0][0], len(w)
    else:
        coords_opt = np.where(s[len(v), :] == score_opt)
        i_opt, j_opt = len(v), coords_opt[0][0]

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
    f = open('rosalind_ba5i.txt', 'r')
    v = f.readline().strip()
    w = f.readline().strip()
    f.close()
    return v, w


v = 'GCACTT'
w = 'CCCAAT'

(score, v_aligned, w_aligned) = overlap_alignment(v, w)
print(score)
print(v_aligned)
print(w_aligned)