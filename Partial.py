from lingpy import *
from lingpy.sequence.sound_classes import *
import networkx as nx

from similarity_networks import local_score

wl = Wordlist('d_chinese.tsv')

def split_morphemes(word):

    out = [[]]
    tokens = [w for w in word]
    while tokens:
        token = tokens.pop(0)
        if token == rc('morpheme_separator'):
            out += [[]]
        else:
            out[-1] += [token]
    return out

mytaxa = [
    "Wenzhou",
    "Guangzhou",
    "Nanning",
    "Xianggang",
    "Haikou",
    "Shantou",
    "Meixian",
    "Taoyuan",
    "Fuzhou",
    "Haerbin",
    "Beijing",
    "Tianjin",
    "Guiyang",
    "Yinchuan",
    "Zhengzhou",
    "Yinchuan",
    "Wulumuqi",
    "Nanjing",
    "Qingdao",
    "Nanchang",
    "Lanzhou",
    "Pingyao"
    ]

idxs = wl.get_list(concept='face', flat=True)
taxa = [wl[idx,'taxon'] for idx in idxs]
words = [wl[idx,'ipa'] for idx in idxs]
chars = [wl[idx,'proto'] for idx in idxs]
morps = [split_morphemes(syllabify(ipa2tokens(w, merge_vowels=False))) for w in words]

def similarity_network(words, threshold):
    """
    Create a similarity network based on partial cognacy scores.
    """
    
    G = nx.Graph()

    m = Multiple(words)
    m.lib_align()


    # attach stuff to alignment
    alm = []


    for i,w1 in enumerate(words):
        if chars[i] in ['面', '脸面', '脸']:
            alm += [[taxa[i],w1]]
            for j,w2 in enumerate(words):
                if i < j:
                    if chars[j] in ['面', '脸面', '脸']:
                        d = local_score(w1, w2) #p.alignments[0][-1]

                        n1 = '{0}'.format(i)
                        n2 = '{0}'.format(j)
                        if n1 not in G:
                            G.add_node(n1,
                                    name=''.join(tokens2class(ipa2tokens(w1), 'sca')),
                                    word=i,
                                    char = chars[i]
                                    )
                        if n2 not in G:
                            G.add_node(n2,
                                name=''.join(tokens2class(ipa2tokens(w2),'sca')), 
                                word=j,
                                char = chars[j]
                                )
                        
                        if d <= threshold:
                            p = Pairwise(w1, w2)
                            p.align(mode='global')
                            almA = tuple(p.alignments[0][0])
                            almB = tuple(p.alignments[0][1])
                            p.align(mode='local')
                            locA = tuple(p.alignments[0][0])
                            locB = tuple(p.alignments[0][1])

                            if locA == almA and locB == almB:
                                eq = 'full'
                            else:
                                eq = 'half'

                            G.add_edge(n1, n2, distance=d, similarity=1-d,
                                    etype=eq)
    
    return G, alm



def partial_cognacies(words, threshold):

    # iterate stuff
    G = nx.Graph()
    
    # make the distance calculation for each segment
    for i,w1 in enumerate(words):
        for j,w2 in enumerate(words):
            if i < j:
                for k,m1 in enumerate(w1):
                    for l,m2 in enumerate(w2):
                        p = Pairwise(m1, m2)
                        p.align(distance=True)
                        d = p.alignments[0][-1]

                        n1 = '{0}_{1}'.format(i,k)
                        n2 = '{0}_{1}'.format(j,l)
                        if n1 not in G:
                            G.add_node(n1,
                                    name=''.join(tokens2class(m1, 'sca')),
                                    word=i,
                                    position=k
                                    )
                        if n2 not in G:
                            G.add_node(n2,
                                name=''.join(tokens2class(m2,'sca')), 
                                word=j,
                                position=l)
                        
                        if d <= threshold:
                            G.add_edge(n1, n2, distance=d)
    
    # check that each segment does only connect to one other segment of another
    # string
    for node,data in G.nodes(data=True):
        
        MN = {}

        # get nodes linking to two sequences
        for neighbor in G.edge[node]:
            
            # get the word
            word = neighbor.split('_')[-1]
            
            try:
                MN[word] += [neighbor]
            except KeyError:
                MN[word] = [neighbor]

        # iterate over keys with more than one neighbor
        for k,v in MN.items():
            if len(v) > 1:
                minD = 1.0

                for neighbor in v:
                    thisD = G.edge[node][neighbor]['distance']
                    if thisD < minD:
                        minD = thisD
                for neighbor in v:
                    thisD = G.edge[node][neighbor]['distance']
                    if thisD > minD:
                        G.remove_edge(node, neighbor)

    return G



def partial_cognacy(w1, w2, threshold, write_gml=False, filename='dummy'):

    m1 = split_morphemes(w1)
    m2 = split_morphemes(w2)

    nodes = []
    for i,m in enumerate(m1):
        nodes += [(i,1)]
    for i,m in enumerate(m2):
        nodes += [(i,2)]
    
    G = nx.Graph()
    for n in nodes:
        G.add_node('{0}_{1}'.format(n[0], n[1]),
                name = ' '.join(tokens2class(w1, 'sca'))
                )

    for i,mA in enumerate(m1):
        for j,mB in enumerate(m2):
            p = Pairwise(mA, mB)
            p.align(distance=True)
            d = p.alignments[0][-1]

            if d <= threshold:
                G.add_edge(
                        '{0}_{1}'.format(i, 1),
                        '{0}_{1}'.format(j, 2),
                        distance = d
                        )
    if write_gml:
        nx.write_gml(G, filename+'.gml')

    # cut the worst nodes 
    for node, data in G.nodes(data=True):

        distance = 1
        for neighbor in G.edge[node]:
            d = G.edge[node][neighbor]['distance']
            if d < distance:
                distance = d
        

        for neighbor in sorted(G.edge[node]):
            if G.edge[node][neighbor]['distance'] > distance:
                G.remove_edge(node, neighbor)

    # now convert stuff to some form
    return G


g, alm = similarity_network(words, 0.3)
nx.write_gml(g, 'r_tmp.gml')

m = Multiple([l[1] for l in alm], merge_vowels=False, model='asjp')
m.prog_align()
with open('r_alignment.msa', 'w') as f:
    f.write('MSA\nFACE\n')
    for i,msa in enumerate(m.alm_matrix):

        f.write(alm[i][0]+'\t'+'\t'.join(msa)+'\n')


import html.parser
test = open('r_tmp.gml').read()

with open('r_face.gml', 'w') as f:
    f.write(html.parser.unescape(test))


