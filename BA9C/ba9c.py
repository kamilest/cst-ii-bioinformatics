def construct_suffix_trie(text):
    trie = []
    trie.append({})

    index = 0
    while len(text) > 0:
        current = 0
        for i in range(len(text)):
            if text[i] in trie[current]:
                current = trie[current][text[i]]
            else:
                # move any other nodes if they start with the same character
                for key in trie[current].keys(): 
                    if key[0] == text[i]:
                        if text[i] not in trie[current]:
                            trie.append({})
                        

                trie[current][text[i:]] = len(trie)
                trie.append({})
                break
        text = text[1:]
        index += 1

    return trie
    


f = open('rosalind_ba9b.txt', 'r')
text = f.readline().strip()
f.close()

suffix_trie = construct_suffix_trie(text)
