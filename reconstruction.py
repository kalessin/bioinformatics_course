from collections import defaultdict, OrderedDict
from operator import itemgetter


def decomposition(k, text):
    """
    decompose a string into a sequence of overlapping k-mers
    >>> decomposition(5, 'CAATCCAAC')
    ['CAATC', 'AATCC', 'ATCCA', 'TCCAA', 'CCAAC']
    >>> decomposition(5, 'AAAAAACGAT')
    ['AAAAA', 'AAAAA', 'AAAAC', 'AAACG', 'AACGA', 'ACGAT']
    """
    n = len(text)
    result = []
    for i in xrange(n - k + 1):
        result.append(text[i:i+k])
    return result


def compose_from_sorted_kmers(kmers):
    """
    compose original string from overlapping k-mers decomposition sequence
    >>> compose_from_sorted_kmers(['ACCGA', 'CCGAA', 'CGAAG', 'GAAGC', 'AAGCT'])
    'ACCGAAGCT'
    """
    result, kmers = kmers[0], kmers[1:]
    for kmer in kmers:
        k = len(kmer)
        assert kmer[0:k-1] == result[-k+1:]
        result += kmer[-1]
    return result


def overlap_graph(kmers):
    """
    For each k-mer in a sequence, generate a list of overlapping k-mers from it
    >>> overlap_graph(['ATGCG', 'GCATG', 'CATGC', 'AGGCA', 'GGCAT'])
    [('GCATG', ['CATGC']), ('CATGC', ['ATGCG']), ('AGGCA', ['GGCAT']), ('GGCAT', ['GCATG'])]
    >>> overlap_graph(['AAAAA', 'AAAAA', 'AAAAC', 'AAACG', 'AACGA', 'ACGAT'])
    [('AAAAA', ['AAAAA', 'AAAAC']), ('AAAAA', ['AAAAA', 'AAAAC']), ('AAAAC', ['AAACG']), ('AAACG', ['AACGA']), ('AACGA', ['ACGAT'])]
    >>> overlap_graph(['AAAAA', 'AAAAA', 'AAAAC', 'AAACG', 'AACGA', 'ACGAA', 'CGAAA', 'GAAAA', 'AAAAA', 'AAAAT'])
    [('AAAAA', ['AAAAA', 'AAAAC', 'AAAAA', 'AAAAT']), ('AAAAA', ['AAAAA', 'AAAAC', 'AAAAA', 'AAAAT']), ('AAAAC', ['AAACG']), ('AAACG', ['AACGA']), ('AACGA', ['ACGAA']), ('ACGAA', ['CGAAA']), ('CGAAA', ['GAAAA']), ('GAAAA', ['AAAAA', 'AAAAA', 'AAAAC', 'AAAAA', 'AAAAT']), ('AAAAA', ['AAAAA', 'AAAAA', 'AAAAC', 'AAAAT'])]
    """
    result = []
    for i in range(len(kmers)):
        kmer = kmers[i]
        k = len(kmer)
        ovp = []
        for j in range(len(kmers)):
            if i != j:
                kkmer = kmers[j]
                if kkmer[0:k-1] == kmer[-k+1:]:
                    ovp.append(kkmer)
        if ovp:
            result.append((kmer, ovp))
    return result


def debruijn_graph(k, text):
    """
    Builds the de Bruijn graph of a sequence
    >>> debruijn_graph(5, 'AAAAAACGAT')
    [('AAAA', ['AAAA', 'AAAA', 'AAAC']), ('AAAC', ['AACG']), ('AACG', ['ACGA']), ('ACGA', ['CGAT'])]
    >>> debruijn_graph(5, 'AGAAAACGAT')
    [('AAAA', ['AAAC']), ('AAAC', ['AACG']), ('AACG', ['ACGA']), ('ACGA', ['CGAT']), ('AGAA', ['GAAA']), ('GAAA', ['AAAA'])]
    >>> debruijn_graph(4, 'AAGATTCTCTAAGA')
    [('AAG', ['AGA', 'AGA']), ('AGA', ['GAT']), ('ATT', ['TTC']), ('CTA', ['TAA']), ('CTC', ['TCT']), ('GAT', ['ATT']), ('TAA', ['AAG']), ('TCT', ['CTA', 'CTC']), ('TTC', ['TCT'])]
    >>> debruijn_graph(2, 'TAATGCCATGGGATGTT')
    [('A', ['A', 'T', 'T', 'T']), ('C', ['A', 'C']), ('G', ['A', 'C', 'G', 'G', 'T']), ('T', ['A', 'G', 'G', 'G', 'T'])]
    """
    kmers = decomposition(k, text)
    return debruijn_graph_from_kmers(k, kmers)


def debruijn_graph_from_kmers(k, kmers):
    """
    de Bruijn graph from ramdonly ordered kmers
    >>> import random
    >>> results = []
    >>> for k, text in [(5, 'AAAAAACGAT'), (5, 'AGAAAACGAT'), (4, 'AAGATTCTCTAAGA'), (2, 'TAATGCCATGGGATGTT')]:
    ...     kmers = decomposition(k, text)
    ...     random.shuffle(kmers)
    ...     results.append(debruijn_graph_from_kmers(k, kmers) == debruijn_graph(k, text))
    >>> all(results)
    True
    """
    debruijn_nodes = defaultdict(list)

    for kmer in kmers:
        debruijn_nodes[kmer[:k-1]].append(kmer[1:])

    return [(d, sorted(dd)) for d, dd in sorted(debruijn_nodes.items(), key=itemgetter(0))]


def walk_cycle(adjacency_map, starting_node=None):
    """
    Find the first cycle in the given eulerian adjacency list
    >>> walk_cycle(OrderedDict([(0, [3]), (1, [0]), (2, [1, 6]), (3, [2]), (4, [2]), (5, [4]), (6, [5, 8]), (7, [9]), (8, [7]), (9, [6])]))
    [0, 3, 2, 6, 8, 7, 9, 6, 5, 4, 2, 1, 0]
    >>> walk_cycle(OrderedDict([(0, [3]), (1, [0]), (2, [6, 1]), (3, [2]), (4, [2]), (5, [4]), (6, [8, 5]), (7, [9]), (8, [7]), (9, [6])]))
    [0, 3, 2, 1, 0]
    """
    if starting_node is not None:
        result = [starting_node]
    else:
        for node in adjacency_map.iterkeys():
            result = [node]
            break
    while True:
        next_nodes = adjacency_map[result[-1]]
        result.append(next_nodes.pop())
        if not next_nodes:
            adjacency_map.pop(result[-2])
        if result[-1] == result[0]:
            break
    return result


def euler_path(adjacency_map):
    """
    >>> euler_path(OrderedDict([(0, [3]), (1, [0]), (2, [6, 1]), (3, [2]), (4, [2]), (5, [4]), (6, [8, 5]), (7, [9]), (8, [7]), (9, [6])]))
    [6, 5, 4, 2, 1, 0, 3, 2, 6, 8, 7, 9, 6]
    """
    cycle = walk_cycle(adjacency_map)
    while True:
        for i, node in enumerate(cycle):
            if node in adjacency_map:
                cycle = cycle[i:-1] + cycle[0:i]
                cycle.extend(walk_cycle(adjacency_map, starting_node=node))
                break
        else:
            break
    return cycle
