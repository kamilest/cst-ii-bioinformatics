def construct_trie(patterns):
    trie = []
    trie.append({})

    for pattern in patterns:
        current = 0
        for i in range(len(pattern)):
            if pattern[i] in trie[current]:
                current = trie[current][pattern[i]]
            else:
                trie[current][pattern[i]] = len(trie)
                current = len(trie)
                trie.append({})

    return trie

def print_trie(trie):
    for i, neighbours in enumerate(trie):
        for c in neighbours:
            print('{}->{}:{}'.format(i, neighbours[c], c))


f = open("rosalind_ba9a.txt", "r")

patterns = []
for line in f:
    patterns.append(line.strip())
f.close()

trie = construct_trie(patterns)
print_trie(trie)
