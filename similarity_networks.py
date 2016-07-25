import networkx as nx
from lingpy import *

def similarity_network(taxa, words, groups, modes='global,local,overlap',
        threshold=0.5, output=False):

    G = nx.Graph()
    for i,(w1, t1, g1) in enumerate(zip(words, taxa, groups)):
        
        cls1 = ''.join(tokens2class(ipa2tokens(w1), 'sca'))
        G.add_node(
                'node_'+str(i+1),
                taxon = t1,
                group = g1,
                name = cls1
                )

        for j,(w2, t2, g2) in enumerate(zip(words, taxa, groups)):
            if i < j:
                distances = {}
                for m in modes.split(','):
                    pair = Pairwise(w1, w2)
                    pair.align(distance=True, mode=m)
                    dst = pair.alignments[0][-1]
                    if dst < 0.5:
                        distances[m] = dst
                if distances:
                    G.add_edge('node_'+str(i+1), 'node_'+str(j+1), **distances)
    
    if output:
        nx.write_gml(G, output+'.gml')
        
    return G

def local_score(a,b):

    p = Pairwise(a,b)
    p.align(mode='local')
    #print(p.alignments, p._alignments, p.tokens[0])
    
    simAB = p.alignments[0][-1]

    seqA = p.alignments[0][0]
    seqB = p.alignments[0][1]
    
    p = Pairwise(seqA,seqA)
    p.align()
    simA = p.alignments[0][-1]

    p = Pairwise(seqB,seqB)
    p.align()
    simB = p.alignments[0][-1]

    # check for a vowel in the alignment
    classes = tokens2class(seqA, 'cv')
    if 'V' in classes and len(seqA) >= 2:

        return 1 - ((2 * simAB)/(simA + simB))
    
    else:
        return 1
    

def normal_score(a,b):

    p = Pairwise(a,b)
    p.align(mode='local', distance=True)
    return p.alignments[0][-1]
