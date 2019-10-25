import numpy as np

def longest_common_subsequence(v, w):
    v_length = len(v)
    w_length = len(w)

    s = np.zeros((v_length+1, w_length+1))
    backtrack = {}
    
    # Dynamic programming
    for i in range(1, v_length+1):
        for j in range(1, w_length+1):
            # Going down or right
            s[i, j] = max(s[i-1, j], s[i, j-1])

            # If match found
            if v[i-1] == w[j-1]:
                s[i, j] = max(s[i, j], s[i-1, j-1]+1)
            
            # Backtrack to previous coordinate
            if s[i, j] == s[i-1, j]:
                backtrack[(i, j)] = (i-1, j)
            elif s[i, j] == s[i, j-1]:
                backtrack[(i, j)] = (i, j-1)
            elif s[i, j] == s[i-1, j-1]+1 and v[i-1] == w[j-1]:
                backtrack[(i, j)] = (i-1, j-1)

    return backtrack

def get_lcs(backtrack, v, w):
    i, j = len(v), len(w)
    out = ''
    while i > 0 and j > 0:
        if backtrack[(i, j)] == (i-1, j):
            i-=1
        elif backtrack[(i, j)] == (i, j-1):
            j-=1
        else:
            i, j = i-1, j-1
            out = v[i] + out

    return out

f = open("rosalind_ba5c.txt", "r")
v = f.readline()
w = f.readline()

backtrack = longest_common_subsequence(v, w)
lcs = get_lcs(backtrack, v, w)
print(lcs)