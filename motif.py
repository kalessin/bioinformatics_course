import random

from math import log
import numpy as np

from frequentwords import computingfrequenciesII, neighbors, numbertopattern
from randomW import randomW


_SYMBOLS = 'ACGT'


def profileMatrix(strings):
    t = len(strings)
    k = len(strings[0])
    count = np.zeros((4, k))
    for string in strings:
        for idx, symbol in enumerate(string):
            count[_SYMBOLS.find(symbol), idx] += 1
    return count / t


def score(strings, entropy=False):
    """
    Returns score of a set of k-mers motifs, and the profile matrix. If entropy is False, the score is the number of mismatches against the
    consensus string.
    """
    profile = profileMatrix(strings)
    if entropy:
        f = np.vectorize(lambda x: 0 if np.isnan(x) else x)
        return np.sum(f(- profile * np.log2(profile))), profile
    else:
        t = len(strings)
        return t * (np.sum(profile) - np.sum(np.max(profile, 0))), profile


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
    that minimizes the hamming distance of pattern to them. It also returns the resulting sum of hamming distances.
    >>> motifspattern('AAA', ['AAAC', 'AAAG', 'AAAT', 'AAAA'])
    (['AAA', 'AAA', 'AAA', 'AAA'], 0)
    >>> motifspattern('AAA', ['AAGC', 'CAAG', 'AAAT', 'AAAA'])
    (['AAG', 'CAA', 'AAA', 'AAA'], 2)
    >>> motifspattern('AAA', ['AAGC', 'CAAG', 'CCGA', 'CGTT'])
    (['AAG', 'CAA', 'CGA', 'CGT'], 7)
    >>> motifspattern('AAA', ['TTACCTTAAC', 'GATATCTGTC', 'ACGGCGTTCG', 'CCCTAAAGAG', 'CGTCAGAGGT'])
    (['TAA', 'ATA', 'ACG', 'AAA', 'AGA'], 5)
    """
    d = 0
    t = len(strings)
    result_motifs = [pattern if pattern in s else None for s in strings]
    remaining = [strings[i] if result_motifs[i] is None else None for i in range(t)]
    result_distance = 0
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


def medianstring(dna, k):
    """
    Generates all possible patterns of length k and searches all them in the target dna strings. Returns
    the ones that matches in all them with the minimal total hamming distance
    >>> medianstring(["AAATTGACGCAT", "GACGACCACGTT", "CGTCAGCGCCTG", "GCTGAGCACCGG", "AGTTCGGGACAG"], 3)
    (['GAC'], 2)
    """
    distance = float('inf')
    medians = []
    for i in xrange(4 ** k):
        pattern = numbertopattern(i, k)
        d = motifspattern(pattern, dna)[1]
        if d < distance:
            distance = d
            medians = [pattern]
        elif d == distance:
            medians.append(pattern)
    return medians, distance


_SYMBOLS = 'ACGT'
def mostprobablekmer(text, k, profile):
    """
    >>> import numpy as np
    >>> text = "ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT"
    >>> k = 5
    >>> profile = np.array([[0.2, 0.2, 0.3, 0.2, 0.3], [0.4, 0.3, 0.1, 0.5, 0.1], [0.3, 0.3, 0.5, 0.2, 0.4], [0.1, 0.2, 0.1, 0.1, 0.2]])
    >>> mostprobablekmer(text, k, profile)
    ('CCGAG', 0.0048000000000000004)
    """
    pmax = 0
    maxkmer = None
    for i in range(len(text) - k + 1):
        kmer = text[i:i+k]
        p = 1
        for j, s in enumerate(kmer):
            ii = _SYMBOLS.find(s)
            p *= profile[ii, j]
        if p > pmax or maxkmer is None:
            pmax = p
            maxkmer = kmer
    return maxkmer, pmax


def profilerandomkmer(text, k, profile):
    """
    Instead of retrieving the most probable k-mer, it picks one randomly, with a probability distribution equal to the probability of each one
    """
    parray = []
    for i in range(len(text) - k + 1):
        kmer = text[i:i+k]
        p = 1
        for j, s in enumerate(kmer):
            i = _SYMBOLS.find(s)
            p *= profile[i, j]
        parray.append(p)
    i = randomW(np.array(parray))
    return text[i:i+k], parray[i]


def mostprobablekmers(dna, k, profile):
    return [mostprobablekmer(s, k, profile)[0] for s in dna]


def greedymotifsearch(dna, k, entropy=False, succession=1):
    """
    Returns a set of motifs, one motif for each input dna string, that minimizes the score function, with
    some biological heuristics for avoiding to search every possible pattern as medianstring does.
    >>> greedymotifsearch(['GGCGTTCAGGCA', 'AAGAATCAGTCA', 'CAAGGAGTTCGC', 'CACGTCAATCAC', 'CAATAATATTCG'], 3, succession=0)
    (['CAG', 'CAG', 'CAA', 'CAA', 'CAA'], 1.9999999999999996)
    >>> greedymotifsearch(['GGCGTTCAGGCA', 'AAGAATCAGTCA', 'CAAGGAGTTCGC', 'CACGTCAATCAC', 'CAATAATATTCG'], 3)
    (['TTC', 'ATC', 'TTC', 'ATC', 'TTC'], 1.9999999999999996)
    >>> greedymotifsearch(['GCCCAA', 'GGCCTG', 'AACCTA', 'TTCCTT'], 3, succession=0)
    (['GCC', 'GCC', 'AAC', 'TTC'], 4.0)
    >>> greedymotifsearch(['GCCCAA', 'GGCCTG', 'AACCTA', 'TTCCTT'], 3)
    (['CCA', 'CCT', 'CCT', 'CCT'], 1.0)
    """
    t = len(dna)
    bestMotifs = [d[0:k] for d in dna]
    bestScore = score(bestMotifs, entropy)[0]
    for i in range(len(dna[0]) - k + 1):
        motifs = [dna[0][i:i+k]]
        sc, profile = score(motifs, entropy)
        for j in range(1, t):
            profile += succession # Laplace's succession rule
            motifs.append(mostprobablekmer(dna[j], k, profile)[0])
            sc, profile = score(motifs, entropy)
        if sc < bestScore:
            bestMotifs = motifs
            bestScore = sc
    return bestMotifs, bestScore


def pickbestrandommotifs(dna, k, n=1, entropy=False):
    bestmotifs = None
    bestscore = float('inf')
    bestprofile = None
    for i in range(n):
        motifs = []
        for string in dna:
            i = random.randint(0, len(string) - k)
            motifs.append(string[i:i + k])
        sc, profile = score(motifs, entropy)
        if sc < bestscore:
            bestmotifs = motifs
            bestscore = sc
            bestprofile = profile
    return bestmotifs, bestscore, bestprofile


def randomizedmotifsearch(dna, k, entropy=False, succession=1):
    """
    >>> while True:
    ...    motifs, score = randomizedmotifsearch(['GCCCAA', 'GGCCTG', 'AACCTA', 'TTCCTT'], 3)
    ...    if score == 1.0:
    ...        break
    >>> motifs in (['CCA', 'CCT', 'CCT', 'CCT'], ['CCC', 'CCT', 'CCT', 'CCT'])
    True
    """
    bestmotifs, bestscore, profile = pickbestrandommotifs(dna, k, entropy=entropy)
    t = len(dna)
    while True:
        profile += succession
        motifs = mostprobablekmers(dna, k, profile)
        sc, profile = score(motifs, entropy)
        if sc < bestscore:
            bestmotifs = motifs
            bestscore = sc
        else:
            return bestmotifs, bestscore


def gibbssampler(dna, k, entropy=False, initmotifs_iters=5, succession=1):
    motifs = list(randomizedmotifsearchx(dna, k, initmotifs_iters, entropy=entropy, succession=succession)[0][0])
    profile = score(motifs, entropy)[1]
    # motifs, _, profile = pickbestrandommotifs(dna, k, 20, entropy=entropy)
    t = len(dna)
    i = random.randint(0, t - 1)
    replaced_motif = profilerandomkmer(dna[i], k, profile)[0]
    motifs = motifs[:i] + [replaced_motif] + motifs[i+1:]
    return motifs, score(motifs, entropy)[0]


def randomizedmotifsearchx(dna, k, iterations=1000, use_gibbssampler=False, entropy=False, succession=1):
    """
    >>> randomizedmotifsearchx(['GCCCAA', 'GGCCTG', 'AACCTA', 'TTCCTT'], 3)
    ([('CCA', 'CCT', 'CCT', 'CCT'), ('CCC', 'CCT', 'CCT', 'CCT')], 1.0)
    >>> while True:
    ...     best, score = randomizedmotifsearchx(['CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA', 'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG', 'TAGTACCGAGACCGAAAGAAGTATACAGGCGT', 'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC', 'AATCCACCAGCTCCACGTGCAATGTTGGCCTA'], 8)
    ...     if score == 9.0:
    ...         break
    >>> set(best).issubset({('AACGGCCA', 'AAGTGCCA', 'TAGTACCG', 'AAGTTTCA', 'ACGTGCAA'), ('TCTCGGGG', 'CCAAGGTG', 'TACAGGCG', 'TTCAGGTG', 'TCCACGTG')})
    True
    """
    best_score = float('inf')
    all_best_motifs = set()
    rmotifsearch = gibbssampler if use_gibbssampler else randomizedmotifsearch
    while iterations > 0:
        iterations -= 1
        motifs, sc = rmotifsearch(dna, k, entropy=entropy, succession=succession)
        if sc < best_score:
            all_best_motifs = {tuple(motifs)}
            best_score = sc
        elif sc == best_score:
            all_best_motifs.add(tuple(motifs))
    return sorted(all_best_motifs), best_score
