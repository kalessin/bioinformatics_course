from collections import defaultdict
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
    debruijn_nodes = defaultdict(list)

    for kmer in kmers:
        debruijn_nodes[kmer[:k-1]].append(kmer[1:])

    return [(d, sorted(dd)) for d, dd in sorted(debruijn_nodes.items(), key=itemgetter(0))]
