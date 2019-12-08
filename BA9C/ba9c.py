def construct_suffix_trie(text):
    # list indexed by node label
    # node label -> (edge label -> (position * next node label))
    trie = []
    trie.append({})

    for i in range(len(text)):
        current = 0
        for j in range(i, len(text)):
            # if there is an outgoing edge from currentNode labeled by currentSymbol
            if text[j] in trie[current]:
                current = trie[current][text[j]][0]
            else:
                trie[current][text[j]] = (len(trie), j, 1)
                trie.append({})
                break
        # if currentNode is a leaf in Trie
        if not trie[current]:
            trie[current][''] = (i, None, None)

    
    # SUFFIXTREECONSTRUCTION(Text)
    # Trie MODIFIEDSUFFIXTRIECONSTRUCTION(Text) for each non-branching path Path in Trie
    # substitute Path by a single edge e connecting the first and last nodes of Path 
    # POSITION(e) <- POSITION(first edge of Path)
    # LENGTH(e) <- number of edges of Path
    # return Trie

    # collapse edges where nodes don't branch
    for node in trie:
        current_node = node
        while len(trie[current_node].keys()) == 1 and \
            '' not in trie[current_node]:
            key = trie[current_node].keys()[0]
            trie[current_node][key][1] = trie[node][0]
            trie[current_node][key][2] = trie[current_node]
            current_node = trie[key][1]



    return trie
    


f = open('ba9c.txt', 'r')
text = f.readline().strip()
f.close()

suffix_trie = construct_suffix_trie(text)
