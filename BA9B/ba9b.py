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


def match_prefix_trie(text, trie):
    i = 0
    v = 0
    while True:
        # if v is a leaf in trie
        if not trie[v]:
            return True
        elif i < len(text) and text[i] in trie[v]:
            v = trie[v][text[i]]
            i+=1
        else:
            return False


def match_trie(text, trie):
    i = 0
    positions = []
    while len(text) > 0:
        if match_prefix_trie(text, trie):
            positions.append(i)
        text = text[1:]
        i+=1
    return positions


f = open('rosalind_ba9b.txt', 'r')
text = f.readline().strip()
patterns = []
for line in f:
    patterns.append(line.strip())
f.close()

trie = construct_trie(patterns)
positions = match_trie(text, trie)
print(*positions)
