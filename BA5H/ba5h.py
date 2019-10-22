import numpy as np

def fitting_alignment(v, w):
    s = np.zeros((len(v)+1, len(w)+1))
    backtrack = {}

    # Initialise scores and backtrack
    for i in range(1, len(v)+1):
        # Scores are already 0.
        backtrack[(i, 0)] = (0, 0)
        
    for j in range(1, len(w)+1):
        s[0, j] = s[0, j-1] - 1
        backtrack[(0, j)] = (0, j-1)

    # Dynamic programming    
    for i in range(1, len(v)+1):
        for j in range(1, len(w)+1):
            # Match or mismatch
            score = 1 if v[i-1] == w[j-1] else -1

            # Indels
            s[i, j] = np.amax([s[i-1, j] - 1, \
                               s[i, j-1] - 1, \
                               s[i-1, j-1] + score])
            
    
            # Backtrack to previous coordinate
            if s[i, j] == s[i-1, j-1] + score:
                backtrack[(i, j)] = (i-1, j-1)
            elif s[i, j] == s[i-1, j] - 1:
                backtrack[(i, j)] = (i-1, j)
            elif s[i, j] == s[i, j-1] - 1:
                backtrack[(i, j)] = (i, j-1)

                
    # Find optimal score and one of the optimal coords.
    score_opt = np.amax(s[:, len(w)])
    coords_opt = np.where(s[:, len(w)] == score_opt)
    (i_opt, j_opt) = coords_opt[0][0], len(w)

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

# Parse the input
f = open('rosalind_ba5h.txt', 'r')
v = f.readline().strip()
w = f.readline().strip()
f.close()

(score, v_aligned, w_aligned) = fitting_alignment(v, w)
print(score)
print(v_aligned)
print(w_aligned)