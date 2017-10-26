"""
Returns the most frequent k-mer patterns in a text, along with the count of ocurrencies
"""
import argparse
from itertools import product
from collections import defaultdict


import patterncount
import hamming


def frequentwords(text, k):
    frequentpatterns = set()
    count = []
    n = len(text)
    for i in range(0, n - k + 1):
        pattern = text[i:i+k]
        count.append(patterncount.patterncount(text, pattern))
    maxcount = max(count)
    if maxcount == 1:
        return set(), 1
    for i in range(0, n - k + 1):
        if count[i] == maxcount:
            frequentpatterns.add(text[i:i+k])
    return frequentpatterns, maxcount


_SYMBOLS = 'ACGT'
def patterntonumber(pattern):
    count = 0
    for c in product(*(_SYMBOLS,) * len(pattern)):
        if ''.join(c) == pattern:
            return count
        count += 1
    raise ValueError(str(pattern))


def patterntonumberR(pattern):
    if not pattern:
        return 0
    symbol = pattern[-1]
    prefix = pattern[:-1]
    return  4 * patterntonumberR(prefix) + _SYMBOLS.find(symbol)


def patterntonumberB(pattern):
    return int(''.join([str(_SYMBOLS.find(c)) for c in pattern]), )


def numbertopattern(i, k):
    count = 0
    for c in product(*(_SYMBOLS,) * k):
        if count == i:
            return ''.join(c)
        count += 1


def numbertopatternR(i, k):
    pattern = ''
    while i > 0:
        i, s = divmod(i, 4)
        pattern = _SYMBOLS[s] + pattern
    while len(pattern) < k:
        pattern = 'A' + pattern
    return pattern


def computingfrequencies(text, k):
    frequencyarray = [0] * 4 ** k
    n = len(text)
    for i in range(0, n - k + 1):
        pattern = text[i:i+k]
        j = patterntonumber(pattern)
        frequencyarray[j] += 1
    return frequencyarray


def fasterfrequentwords(text, k):
    frequentpatterns = set()
    frequencyarray = computingfrequencies(text, k)
    maxcount = max(frequencyarray)
    for i in range(0, 4 ** k - 1):
        if frequencyarray[i] == maxcount:
            frequentpatterns.add(numbertopattern(i, k))
    return frequentpatterns, maxcount


def neighbors(pattern, d):
    """
    >>> neighbors('AA', 1)
    set(['AA', 'AC', 'AG', 'CA', 'AT', 'GA', 'TA'])
    """
    if d == 0:
        return {pattern}
    if len(pattern) == 1:
        return set(_SYMBOLS)
    first, suffix = pattern[0], pattern[1:]
    neighborhood = set()
    for sn in neighbors(suffix, d):
        if hamming.hamming(suffix, sn) < d:
            for symbol in _SYMBOLS:
                neighborhood.add(symbol + sn)
        else: # hamming.hamming(suffix, sn) == d
            neighborhood.add(first + sn)
    return neighborhood


def computingfrequenciesII(text, k, d=0):
    """
    >>> dict(computingfrequenciesII('AAAAAAAAAA', 2))
    {'AA': 9}
    >>> dict(computingfrequenciesII('AAAAAAAAAA', 2, 1))
    {'AA': 9, 'AC': 9, 'AT': 9, 'AG': 9, 'CA': 9, 'GA': 9, 'TA': 9}
    >>> dict(computingfrequenciesII('TAGCG', 2, 1))
    {'AA': 2, 'AC': 2, 'GT': 1, 'AG': 2, 'CC': 2, 'CA': 2, 'CG': 2, 'GG': 3, 'GC': 1, 'AT': 1, 'GA': 2, 'CT': 1, 'TT': 1, 'TG': 3, 'TC': 2, 'TA': 1}
    """
    frequencyarray = defaultdict(int)
    n = len(text)
    for i in range(0, n - k + 1):
        pattern = text[i:i+k]
        frequencyarray[pattern] += 1
    if d > 0: # avoid unneeded step if d == 0
        mfrequencyarray = frequencyarray.copy()
        patterns = frequencyarray.keys()
        for i in range(0, len(patterns)):
            for neighbor in neighbors(patterns[i], d):
                if neighbor != patterns[i]:
                    mfrequencyarray[patterns[i]] += frequencyarray[neighbor]
                    if neighbor not in patterns: # we will count it later
                        mfrequencyarray[neighbor] += frequencyarray[patterns[i]]
        return mfrequencyarray
    return frequencyarray


def frequentwordsII(text, k, d=0):
    """
    The fastest one
    >>> s = frequentwordsII('AGTCAGTC', 4, 2)[0]
    >>> s == {'TCTC', 'CGGC', 'AAGC', 'TGTG', 'GGCC', 'AGGT', 'ATCC', 'ACTG', 'ACAC', 'AGAG', 'ATTA', 'TGAC', 'AATT', 'CGTT', 'GTTC', 'GGTA', 'AGCA', 'CATC'}
    True
    >>> frequentwordsII('AATTAATTGGTAGGTAGGTA', 4)[0]
    set(['GGTA'])
    >>> s = frequentwordsII('ATA', 3, 1)[0]
    >>> s == {'AAA', 'ACA', 'AGA', 'ATA', 'ATC', 'ATG', 'ATT', 'CTA', 'GTA', 'TTA'}
    True
    >>> frequentwordsII('AAT', 3, 0)[0]
    set(['AAT'])
    >>> frequentwordsII('TAGCG', 2, 1)[0]
    set(['GG', 'TG'])
    >>> s = frequentwordsII('ACGTTGCATGTCGCATGATGCATGAGAGCT', 4, 1)[0]
    >>> s == set(['ATGC', 'ATGT', 'GATG'])
    True
    """
    frequencyarray = computingfrequenciesII(text, k, d)
    maxcount = max(frequencyarray.itervalues())
    frequentpatterns = set()
    for p, c in frequencyarray.iteritems():
        if c == maxcount:
            frequentpatterns.add(p)
    return frequentpatterns, maxcount

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("text", metavar="TEXTFILE", type=argparse.FileType("r"))
    parser.add_argument("k", type=int)
    parser.add_argument("--alg", type=int, default=2)
    parser.add_argument("--hamming", type=int, default=0)
    args = parser.parse_args()
    algmap = {
        0: frequentwords,
        1: fasterfrequentwords,
        2: frequentwordsII,
    }

    print "frequent words: %s (%d)" % algmap[args.alg](args.text.read().strip(), args.k, args.hamming)
