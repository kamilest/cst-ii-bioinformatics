import numpy as np

def nussinov_folding(v):
  n = len(v)
  s = np.zeros((n, n), dtype=np.int)

  def score(a, b):
    if a == 'A':
      return b == 'U'
    elif a == 'U':
      return b == 'A'
    elif a == 'C':
      return b == 'G'
    else: # a == 'G'
      return b == 'C'
  
  for i in range(n-1):
    s[i,i] = 0
    s[i+1,i] = 0

  s[n-1,n-1] = 0
  for k in range(1,n):
    for i in range(n-k):
      j = i + k

      if range(i+1, j-1):
        s_bi = [s[i,l] + s[l+1,j] for l in range(i+1, j-1)]
      else:
        s_bi = [0]

      s[i,j] = np.amax([s[i+1,j],\
                        s[i,j-1],\
                        s[i+1,j-1] + score(v[i],v[j]), \
                        np.amax(s_bi)])
  
  print(v)
  print(s)


# Parse the input
f = open('in.txt', 'r')
v = f.readline().strip()
f.close()

nussinov_folding(v)