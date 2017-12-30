from copy import deepcopy
from collections import defaultdict, OrderedDict
from operator import itemgetter
from itertools import permutations


def composition(k, text):
    """
    decompose a string into a sequence of overlapping k-mers
    >>> list(composition(5, 'CAATCCAAC'))
    ['CAATC', 'AATCC', 'ATCCA', 'TCCAA', 'CCAAC']
    >>> list(composition(5, 'AAAAAACGAT'))
    ['AAAAA', 'AAAAA', 'AAAAC', 'AAACG', 'AACGA', 'ACGAT']
    """
    n = len(text)
    for i in xrange(n - k + 1):
        yield text[i:i+k]


def compose_from_sorted_kmers(kmers):
    """
    compose original string from overlapping k-mers composition sequence
    >>> compose_from_sorted_kmers(['ACCGA', 'CCGAA', 'CGAAG', 'GAAGC', 'AAGCT'])
    'ACCGAAGCT'
    """
    result, kmers = kmers[0], kmers[1:]
    k = len(result)
    for kmer in kmers:
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
    kmers = composition(k, text)
    return debruijn_graph_from_kmers(k, kmers)


def debruijn_graph_from_kmers(k, kmers):
    """
    de Bruijn graph from ramdonly ordered kmers
    >>> import random
    >>> results = []
    >>> for k, text in [(5, 'AAAAAACGAT'), (5, 'AGAAAACGAT'), (4, 'AAGATTCTCTAAGA'), (2, 'TAATGCCATGGGATGTT')]:
    ...     kmers = list(composition(k, text))
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


class CloseLink(object):
    pass


close_link = CloseLink()


def euler_path(adjacency_map, path_permutations=0):
    """
    >>> euler_path(OrderedDict([(0, [3]), (1, [0]), (2, [6, 1]), (3, [2]), (4, [2]), (5, [4]), (6, [8, 5]), (7, [9]), (8, [7]), (9, [6])]))
    [6, 5, 4, 2, 1, 0, 3, 2, 6, 8, 7, 9, 6]
    >>> euler_path(OrderedDict([(0, [2]), (1, [3]), (2, [1]), (3, [0, 4]), (6, [3, 7]), (7, [8]), (8, [9]), (9, [6])]))
    [6, 7, 8, 9, 6, 3, 0, 2, 1, 3, 4]
    >>> euler_path(OrderedDict([(0, [2, 3]), (1, [3]), (2, [1]), (3, [0, 4, 5]), (5, [0]), (6, [3, 7]), (7, [8]), (8, [9]), (9, [6])]))
    [6, 7, 8, 9, 6, 3, 0, 3, 5, 0, 2, 1, 3, 4]
    >>> euler_path(OrderedDict([(0, [2, 3]), (1, [3]), (2, [1, 4]), (3, [0, 4, 5]), (4, [2]), (5, [0]), (6, [3, 7]), (7, [8]), (8, [9]), (9, [6])]))
    [6, 7, 8, 9, 6, 3, 4, 2, 1, 3, 0, 3, 5, 0, 2, 4]
    >>> euler_path(OrderedDict([(0, [2, 3]), (1, [3]), (2, [1, 4]), (3, [0, 4, 5]), (4, [2]), (5, [0]), (6, [3, 7]), (7, [8]), (8, [9]), (9, [6])]), 1)
    [6, 7, 8, 9, 6, 3, 5, 0, 3, 4, 2, 1, 3, 0, 2, 4]
    >>> euler_path(OrderedDict([(0, [2, 3]), (1, [3]), (2, [1, 4]), (3, [0, 4, 5]), (4, [2]), (5, [0]), (6, [3, 7]), (7, [8]), (8, [9]), (9, [6])]), 2)
    [6, 7, 8, 9, 6, 3, 0, 2, 4, 2, 1, 3, 5, 0, 3, 4]
    >>> euler_path(OrderedDict([(0, [2, 3]), (1, [3]), (2, [1, 4]), (3, [0, 4, 5]), (4, [2]), (5, [0]), (6, [3, 7]), (7, [8]), (8, [9]), (9, [6])]), 3)
    [6, 7, 8, 9, 6, 3, 5, 0, 3, 0, 2, 4, 2, 1, 3, 4]
    >>> euler_path(OrderedDict([(0, [2, 3]), (1, [3]), (2, [1, 4]), (3, [0, 4, 5]), (4, [2]), (5, [0]), (6, [3, 7]), (7, [8]), (8, [9]), (9, [6])]), 9)
    Traceback (most recent call last):
        ...
    AssertionError
    """
    adjacency_map = deepcopy(adjacency_map)
    for node, next_nodes in adjacency_map.iteritems():
        if path_permutations > 0:
            permuted_next_nodes = permutations(next_nodes)
            permuted_next_nodes.next()
            for next_nodes in permuted_next_nodes:
                path_permutations -= 1
                adjacency_map[node] = list(next_nodes)
                if path_permutations == 0:
                    break
    assert path_permutations == 0
    inverse_adjacency_map = defaultdict(list)
    for node, next_nodes in adjacency_map.iteritems():
        for next_node in next_nodes:
            inverse_adjacency_map[next_node].append(node)
    first_node = last_node = None
    for node, next_nodes in adjacency_map.iteritems():
        if len(next_nodes) > len(inverse_adjacency_map.get(node, [])):
            first_node = node
            break
    for node, next_nodes in inverse_adjacency_map.iteritems():
        if len(next_nodes) > len(adjacency_map.get(node, [])):
            last_node = node
            adjacency_map.setdefault(last_node, []).append(close_link)
            adjacency_map[close_link] = [first_node]
            break
    cycle = walk_cycle(adjacency_map)
    while True:
        for i, node in enumerate(cycle):
            if node in adjacency_map:
                cycle = cycle[i:-1] + cycle[0:i]
                cycle.extend(walk_cycle(adjacency_map, starting_node=node))
                break
        else:
            break
    for i, node in enumerate(cycle):
        if node == close_link:
            cycle = cycle[i+1:] + cycle[1:i]
            break
    return cycle


def sequence_from_kmers(k, kmers):
    """
    >>> sequence_from_kmers(4, ['CTTA', 'ACCA', 'TACC', 'GGCT', 'GCTT', 'TTAC'])
    'GGCTTACCA'
    """
    adjacency_map = OrderedDict(debruijn_graph_from_kmers(k, kmers))
    return compose_from_sorted_kmers(euler_path(adjacency_map))


def universal_string(k):
    kmers = [bin(i)[2:].zfill(k) for i in range(2**k)]
    return sequence_from_kmers(k, kmers)


def universal_circular_string(k):
    kmers = [bin(i)[2:].zfill(k) for i in range(2**k)]
    debruijn = dict(debruijn_graph_from_kmers(k, kmers))
    return compose_from_sorted_kmers(euler_path(debruijn)[:-k+1])


def gapped_pairs_composition(k, d, text):
    """
    >>> list(gapped_pairs_composition(3, 1, 'TAATGCCATGGGATGTT'))
    [('TAA', 'GCC'), ('AAT', 'CCA'), ('ATG', 'CAT'), ('TGC', 'ATG'), ('GCC', 'TGG'), ('CCA', 'GGG'), ('CAT', 'GGA'), ('ATG', 'GAT'), ('TGG', 'ATG'), ('GGG', 'TGT'), ('GGA', 'GTT')]
    """
    queue = []
    nlen = d + k + 1
    for kmer in composition(k, text):
        queue.append(kmer)
        if len(queue) == nlen:
            yield queue.pop(0), kmer


def compose_from_sorted_gapped_pairs(d, paired_kmers):
    """
    >>> compose_from_sorted_gapped_pairs(1, [('TAA', 'GCC'), ('AAT', 'CCA'), ('ATG', 'CAT'), ('TGC', 'ATG'), ('GCC', 'TGG'), ('CCA', 'GGG'), ('CAT', 'GGA'), ('ATG', 'GAT'), ('TGG', 'ATG'), ('GGG', 'TGT'), ('GGA', 'GTT')])
    'TAATGCCATGGGATGTT'
    >>> compose_from_sorted_gapped_pairs(2, [('GACC', 'GCGC'), ('ACCG', 'CGCC'), ('CCGA', 'GCCG'), ('CGAG', 'CCGG'), ('GAGC', 'CGGA')])
    'GACCGAGCGCCGGA'
    >>> compose_from_sorted_gapped_pairs(3, [('GTG', 'GTG'), ('TGG', 'TGA'), ('GGT', 'GAG'), ('GTC', 'AGA'), ('TCG', 'GAT'), ('CGT', 'ATG'), ('GTG', 'TGT'), ('TGA', 'GTT'), ('GAG', 'TTG'), ('AGA', 'TGA')])
    'GTGGTCGTGAGATGTTGA'
    """
    k = len(paired_kmers[0][0])
    resultP = compose_from_sorted_kmers(map(itemgetter(0), paired_kmers))
    resultS = compose_from_sorted_kmers(map(itemgetter(1), paired_kmers))
    if resultP[k+d:] == resultS[:-k-d]:
        return resultP + resultS[-k-d:]


def debruijn_graph_from_gapped_pairs(k, paired_kmers):
    return debruijn_graph_from_kmers(k, map(itemgetter(0), paired_kmers)), debruijn_graph_from_kmers(k, map(itemgetter(1), paired_kmers))


def euler_path_for_gapped_pairs(k, d, prefix_adjacency_map, suffix_adjacency_map):
    permutations = 0
    while True:
        peuler_path = euler_path(prefix_adjacency_map, permutations)
        seuler_path = euler_path(suffix_adjacency_map)
        if peuler_path[k+d:] == seuler_path[:-k-d]:
            break
        else:
            permutations += 1
    return zip(peuler_path, seuler_path)


def sequence_from_gapped_pairs(k, d, paired_kmers):
    """
    >>> sequence_from_gapped_pairs(4, 2, [('GAGA', 'TTGA'), ('TCGT', 'GATG'), ('CGTG', 'ATGT'), ('TGGT', 'TGAG'), ('GTGA', 'TGTT'), ('GTGG', 'GTGA'), ('TGAG', 'GTTG'), ('GGTC', 'GAGA'), ('GTCG', 'AGAT')])
    'GTGGTCGTGAGATGTTGA'
    """
    adjacency_maps = debruijn_graph_from_gapped_pairs(k, paired_kmers)
    epath = euler_path_for_gapped_pairs(k, d, OrderedDict(adjacency_maps[0]), OrderedDict(adjacency_maps[1]))
    return compose_from_sorted_gapped_pairs(d + 1, epath)
