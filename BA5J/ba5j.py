import numpy as np

def affine_gap_alignment(v, w, scoring_matrix, gap_opening_penalty, gap_extension_penalty):
    # Diagonal edges of weight Score(vi, wj) representing matches and mismatches
    middle = np.matrix(np.ones((len(v)+1, len(w)+1)) * -np.inf)
    middle[0, 0] = 0

    # Only horizontal edges with weights -e to represent gap extensions in w
    upper = np.matrix(np.ones((len(v)+1, len(w)+1)) * -np.inf)

    # Only vertical edges with weight -e to represent gap extensions in v
    lower = np.matrix(np.ones((len(v)+1, len(w)+1)) * -np.inf)

    backtrack = {}

    # Initialise first row and column
    for i in range(1, len(v)+1):
        # middle[i, 0] = - gap_opening_penalty - (i-1) * gap_extension_penalty
        lower[i, 0] = - gap_opening_penalty - (i-1) * gap_extension_penalty
        backtrack[(i, 0)] = (i-1, 0)
    
    for j in range(1, len(w)+1):
        # middle[0, j] = - gap_opening_penalty - (j-1) * gap_extension_penalty
        upper[0, j] = - gap_opening_penalty - (j-1) * gap_extension_penalty
        backtrack[(0, j)] = (0, j-1)

    # Dynamic programming
    for i in range(1, len(v)+1):
        for j in range(1, len(w)+1):
            lower[i, j] = max(middle[i-1, j] - gap_opening_penalty, \
                lower[i-1, j] - gap_extension_penalty)
            
            upper[i, j] = max(middle[i, j-1] - gap_opening_penalty, \
                upper[i, j-1] - gap_extension_penalty)

            # If match or mismatch
            score = scoring_matrix[(v[i-1], w[j-1])]
            middle[i, j] = np.amax([middle[i-1, j-1] + score, \
                lower[i, j], upper[i, j]])

            # Backtrack to previous coordinate
            if middle[i, j] == middle[i-1, j-1] + score:
                backtrack[(i, j)] = (i-1, j-1)
            elif middle[i, j] == lower[i, j]:
                backtrack[(i, j)] = (i-1, j)
            elif middle[i, j] == upper[i, j]:
                backtrack[(i, j)] = (i, j-1)
                
    score = middle[len(v), len(w)]
    print(middle)
    print(lower)
    print(upper)

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
    
    return (int(score), v_aligned, w_aligned)


# Parse the input
f = open('ba5j.txt', 'r')
v = f.readline().strip()
w = f.readline().strip()
f.close()

# Parse the scoring matrix
f = open('blosum62.txt', 'r')
aas = f.readline().split() # amino acid codes
scoring_matrix = {}

for i in range(len(aas)):
    scores_i = f.readline().split()
    for j in range(len(scores_i)):
        scoring_matrix[(aas[i], aas[j])] = int(scores_i[j])
f.close()

(score, v_aligned, w_aligned) = affine_gap_alignment(v, w, scoring_matrix, 11, 1)
print(score)
print(v_aligned)
print(w_aligned)