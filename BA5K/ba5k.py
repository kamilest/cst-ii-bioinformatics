import numpy as np

def from_source_lin(v, w, scoring_matrix, indel_penalty):
    middle_column = int(np.floor(len(w) / 2.0))
    s = np.zeros((len(v)+1, 2)) # linear space

    # Scoring from source to middle column computed in linear space.
    for i in range(1, len(v)+1):
        s[i, 0] = s[i-1, 0] - indel_penalty
    
    for col in range(1, middle_column+1):
        s_col = col % 2
        prev_col = (col + 1) % 2
        s[0, s_col] = s[0, prev_col] - indel_penalty
        for i in range(1, len(v)+1):
            score = scoring_matrix[(v[i-1], w[col-1])]
            s[i, s_col] = max([\
                s[i-1, s_col] - indel_penalty, \
                s[i, prev_col] - indel_penalty, \
                s[i-1, prev_col] + score])

    from_source = s[:, middle_column % 2]
    return from_source

def to_sink_lin(v, w, scoring_matrix, indel_penalty):
    middle_column = int(np.ceil(len(w) / 2.0))
    s = np.zeros((len(v)+1, 2)) # linear space

    v_rev = v[::-1]
    w_rev = w[::-1]

    for i in range(1, len(v)+1):
        s[i, 0] = s[i-1, 0] - indel_penalty
    
    backtrack_row = np.zeros(len(v)+1)
    for col in range(1, middle_column+1):
        s_col = col % 2
        prev_col = (col + 1) % 2
        s[0, s_col] = s[0, prev_col] - indel_penalty
        for i in range(1, len(v)+1):
            score = scoring_matrix[(v_rev[i-1], w_rev[col-1])]
            s[i, s_col] = max([\
                s[i-1, s_col] - indel_penalty, \
                s[i, prev_col] - indel_penalty, \
                s[i-1, prev_col] + score])
            if s[i, s_col] == s[i, prev_col] - indel_penalty:
                backtrack_row[i] = i
            else:
                backtrack_row[i] = i-1

    to_sink = np.flip(s[:, middle_column % 2])
    return (to_sink, np.flip(backtrack_row))

def middle_edge_linspace(v, w, scoring_matrix, indel_penalty):
    from_source = from_source_lin(v, w, scoring_matrix, indel_penalty)
    (to_sink, next_i) = to_sink_lin(v, w, scoring_matrix, indel_penalty)

    lengths = from_source + to_sink
    index = np.argmax(lengths)
    next_index = int(len(v) - next_i[index])

    return((index, len(w) // 2), (next_index, len(w) // 2 + 1))


# Parse the input
f = open('rosalind_ba5k.txt', 'r')
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

((i, j), (k, l)) = middle_edge_linspace(v, w, scoring_matrix, 5)
print("({}, {}) ({}, {})".format(i, j, k, l))