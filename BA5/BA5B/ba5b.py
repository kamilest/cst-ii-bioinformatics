"""
http://rosalind.info/problems/ba5b/

Length of a Longest Path in the Manhattan Tourist Problem

Find the length of a longest path in a rectangular city.

Given: Integers n and m, followed by an n × (m+1) matrix Down and an (n+1) × 
m matrix Right. The two matrices are separated by the "-" symbol.

Return: The length of a longest path from source (0, 0) to sink (n, m) in the 
n × m rectangular grid whose edges are defined by the matrices Down and Right.

Sample Dataset

4 4
1 0 2 4 3
4 6 5 2 1
4 4 5 2 1
5 6 8 5 3
-
3 2 4 0
3 2 4 2
0 7 3 3
3 3 0 2
1 3 2 2

Sample Output
34
"""
import numpy as np

def manhattan_tourist(n, m, down, right):
    s = np.zeros((n+1, m+1))

    # Initialise the column downwards
    for i in range(1, n+1):
        s[i, 0] = s[i-1, 0] + down[i-1, 0]
    
    # Initialise the first row to the right
    for j in range(1, m+1):
        s[0, j] = s[0, j-1] + right[0, j-1]
    
    # Dynamic programming
    for i in range(1, n+1):
        for j in range(1, m+1):
            s[i, j] = max(s[i-1, j] + down[i-1, j], \
                 s[i, j-1] + right[i, j-1])

    return s[n, m]

f = open("rosalind_ba5b.txt", "r")
[n, m] = [int(s) for s in f.readline().split(" ")]

down = np.zeros((n, m+1))
right = np.zeros((n+1, m))

for i in range(n):
    down[i, :] = [int(s) for s in f.readline().split(" ")]

f.readline() # read the dash

for i in range(n+1):
    right[i, :] = [int(s) for s in f.readline().split(" ")]

print(manhattan_tourist(n, m, down, right))