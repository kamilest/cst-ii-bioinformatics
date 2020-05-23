import numpy as np

def dna_scoring(a, b, match_score, mismatch_penalty):
    if a == b:
        return match_score
    else:
        return mismatch_penalty


def affine_gap_alignment(v, w, match_score, mismatch_penalty, gap_opening_penalty, gap_extension_penalty, scoring_matrix=None):
    # Diagonal edges of weight Score(vi, wj) representing matches and mismatches
    middle = np.matrix(np.ones((len(v)+1, len(w)+1)) * -np.inf)
    middle[0, 0] = 0

    # Only horizontal edges with weights -e to represent gap extensions in w
    upper = np.matrix(np.ones((len(v)+1, len(w)+1)) * -np.inf)

    # Only vertical edges with weight -e to represent gap extensions in v
    lower = np.matrix(np.ones((len(v)+1, len(w)+1)) * -np.inf)

    backtrack_lower = {}
    backtrack_middle = {}
    backtrack_upper = {}

    # Initialise first row and column
    for i in range(1, len(v)+1):
        middle[i, 0] = - gap_opening_penalty - (i-1) * gap_extension_penalty
        lower[i, 0] = - gap_opening_penalty - (i-1) * gap_extension_penalty
        backtrack_lower[(i, 0)] = 'lower'
    
    for j in range(1, len(w)+1):
        middle[0, j] = - gap_opening_penalty - (j-1) * gap_extension_penalty
        upper[0, j] = - gap_opening_penalty - (j-1) * gap_extension_penalty
        backtrack_upper[(0, j)] = 'upper'

    # Dynamic programming
    for i in range(1, len(v)+1):
        for j in range(1, len(w)+1):
            lower[i, j] = max(middle[i-1, j] - gap_opening_penalty, \
                lower[i-1, j] - gap_extension_penalty)
            if lower[i, j] == middle[i-1, j] - gap_opening_penalty:
                backtrack_lower[(i, j)] = 'middle'
            else:
                backtrack_lower[(i, j)] = 'lower'
            
            upper[i, j] = max(middle[i, j-1] - gap_opening_penalty, \
                upper[i, j-1] - gap_extension_penalty)
            if upper[i, j] == middle[i, j-1] - gap_opening_penalty:
                backtrack_upper[(i, j)] = 'middle'
            else:
                backtrack_upper[(i, j)] = 'upper'

            # If match or mismatch
            if scoring_matrix is not None:
                score = scoring_matrix[(v[i-1], w[j-1])]
            else: 
                score = dna_scoring(v[i-1], w[j-1], match_score, -1 * mismatch_penalty)

            middle[i, j] = np.amax([middle[i-1, j-1] + score, \
                lower[i, j], upper[i, j]])

            # Backtrack to previous coordinate
            if middle[i, j] == middle[i-1, j-1] + score:
                backtrack_middle[(i, j)] = 'middle'
            elif middle[i, j] == lower[i, j]:
                backtrack_middle[(i, j)] = 'lower'
            elif middle[i, j] == upper[i, j]:
                backtrack_middle[(i, j)] = 'upper'
                
    score = middle[len(v), len(w)]
    backtrack_matrix = 'middle'

    # Backtrack to find alignments
    i, j = len(v), len(w)
    v_aligned, w_aligned = '', ''
    while i > 0 and j > 0:
        if backtrack_matrix == 'middle':
            if backtrack_middle[(i, j)] == 'upper':
                backtrack_matrix = 'upper'
            elif backtrack_middle[(i, j)] == 'lower':
                backtrack_matrix = 'lower'
            else:
                i, j = i-1, j-1
                v_aligned = v[i] + v_aligned
                w_aligned = w[j] + w_aligned
        elif backtrack_matrix == 'lower':
            if backtrack_lower[(i, j)] == 'middle':
                backtrack_matrix = 'middle'
            i-=1
            v_aligned = v[i] + v_aligned
            w_aligned = '-' + w_aligned

        elif backtrack_matrix == 'upper':
            if backtrack_upper[(i, j)] == 'middle':
                backtrack_matrix = 'middle'
            j-=1
            v_aligned = '-' + v_aligned
            w_aligned = w[j] + w_aligned
    
    return (int(score), v_aligned, w_aligned)

def parse_input_strings():
    # Parse the input
    f = open('rosalind_ba5j.txt', 'r')
    v = f.readline().strip()
    w = f.readline().strip()
    f.close()
    return v, w

def parse_scoring_matrix():
    # Parse the scoring matrix
    f = open('blosum62.txt', 'r')
    aas = f.readline().split() # amino acid codes
    scoring_matrix = {}

    for i in range(len(aas)):
        scores_i = f.readline().split()
        for j in range(len(scores_i)):
            scoring_matrix[(aas[i], aas[j])] = int(scores_i[j])
    f.close()
    
    return scoring_matrix


match_score = 1
mismatch_penalty = 1
indel_penalty = 4
gap_opening_penalty = 11
gap_extension_penalty=1

assert(match_score >= 0)
assert(mismatch_penalty >= 0)
assert(indel_penalty >= 0)
assert(gap_opening_penalty >= 0)
assert(gap_extension_penalty >= 0)

v = 'GCACTT'
w = 'CCCAAT'

(score, v_aligned, w_aligned) = affine_gap_alignment(v, w, match_score, mismatch_penalty, gap_opening_penalty, gap_extension_penalty)

print(score)
print(v_aligned)
print(w_aligned)