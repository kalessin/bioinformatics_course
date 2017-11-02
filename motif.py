from math import log
import numpy as np

from frequentwords import computingfrequenciesII, neighbors


_SYMBOLS = 'ACGT'


def score(strings):
    """
    Returns entropy of a set of k-mers
    """
    t = len(strings)
    k = len(strings[0])
    count = np.zeros((4, k))
    for string in strings:
        for idx, symbol in enumerate(string):
            count[_SYMBOLS.find(symbol), idx] += 1
    profile = count / t
    f = np.vectorize(lambda x: 0 if np.isnan(x) else x)
    return np.sum(f(- profile * np.log2(profile)))


def motiffind_basic(strings, k, d):
    """
    Motif finder based on getting intersections of each set of all k-mers (with d mismatches)
    from each input string
    >>> motiffind_basic(["ATTTGGC", "TGCCTTA", "CGGTATC", "GAAAATT"], 3, 1)
    set(['ATT', 'TTT', 'GTT', 'ATA'])
    >>> motiffind_basic(["ACGT", "ACGT", "ACGT"], 3, 0)
    set(['ACG', 'CGT'])
    >>> motiffind_basic(["AAAAA", "AAAAA", "AAAAA"], 3, 1)
    set(['AAG', 'AAA', 'AAC', 'ATA', 'ACA', 'AGA', 'AAT', 'TAA', 'CAA', 'GAA'])
    >>> l = sorted(motiffind_basic(["AAAAA", "AAAAA", "AAAAA"], 3, 3))
    >>> len(l) == 64 and set([len(i) for i in l]) == {3}
    True
    >>> motiffind_basic(["AAAAA", "AAAAA", "AACAA"], 3, 0)
    set([])
    >>> motiffind_basic(["AACAA", "AAAAA", "AAAAA"], 3, 0)
    set([])
    """
    first, strings = strings[0], strings[1:]
    patterns = set(computingfrequenciesII(first, k, d))
    for string in strings:
        patterns.intersection_update(computingfrequenciesII(string, k, d).keys())
    return patterns


def motifspattern(pattern, strings):
    """
    Operates in inverse way than motiffind_basic: given a candidate pattern, return a list of motifs, one from each string,
    that minimizes the hamming distance of pattern to them
    >>> motifspattern('AAA', ['AAAC', 'AAAG', 'AAAT', 'AAAA'])
    (['AAA', 'AAA', 'AAA', 'AAA'], 0)
    >>> motifspattern('AAA', ['AAGC', 'CAAG', 'AAAT', 'AAAA'])
    (['AAG', 'CAA', 'AAA', 'AAA'], 2)
    >>> motifspattern('AAA', ['AAGC', 'CAAG', 'CCGA', 'CGTT'])
    (['AAG', 'CAA', 'CGA', 'CGT'], 7)
    """
    d = 0
    t = len(strings)
    result_motifs = [pattern if pattern in s else None for s in strings]
    remaining = [strings[i] if result_motifs[i] is None else None for i in range(t)]
    result_distance = 0
    remaining = [strings[i] for i, s in enumerate(result_motifs) if s is None]
    k = len(pattern)
    processed_neighbors = {s: {pattern} for s in strings}
    while any(remaining):
        d += 1
        for i, s in enumerate(remaining):
            if s is None:
                continue
            candidate_nn = None
            candidate_index = len(s)
            for n in list(processed_neighbors[s]):
                for nn in neighbors(n, 1):
                    if nn not in processed_neighbors[s]:
                        processed_neighbors[s].add(nn)
                        idx = s.find(nn)
                        if 0 <= idx < candidate_index:
                            candidate_nn = nn
                            candidate_index = idx
            if candidate_nn is not None:
                result_motifs[i] = candidate_nn
                result_distance += d
                remaining = [strings[ii] if result_motifs[ii] is None else None for ii in range(t)]
                processed_neighbors.pop(s)
            if not any(remaining):
                return result_motifs, result_distance
    return result_motifs, result_distance
