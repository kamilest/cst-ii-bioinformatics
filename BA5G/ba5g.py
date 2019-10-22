import numpy as np

def edit_distance(v, w):
    s = np.zeros((len(v)+1, len(w)+1))

    # Initialise first row and column
    for i in range(1, len(v)+1):
        s[i, 0] = s[i-1, 0] + 1
    
    for j in range(1, len(w)+1):
        s[0, j] = s[0, j-1] + 1

    # Dynamic programming
    for i in range(1, len(v)+1):
        for j in range(1, len(w)+1):
            score = 0 if v[i-1] == w[j-1] else 1
            s[i, j] = np.amin([s[i-1, j] + 1, \
                               s[i, j-1] + 1, \
                               s[i-1, j-1] + score])
                
    return int(s[len(v), len(w)])

# Parse the input
f = open('rosalind_ba5g.txt', 'r')
v = f.readline().strip()
w = f.readline().strip()
f.close()

print(edit_distance(v, w))