"""
Returns the most frequent k-mer patterns in a text, along with the count of ocurrencies
"""
import argparse
from itertools import product
from collections import defaultdict


import patterncount
import hamming
import rcomplement
import multifile


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


def computingfrequenciesII(text, k, d=0, rc=False, start=0, end=None):
    """
    >>> dict(computingfrequenciesII('AAAAAAAAAA', 2))
    {'AA': 9}
    >>> dict(computingfrequenciesII('AAAAAAAAAA', 2, 1))
    {'AA': 9, 'AC': 9, 'AT': 9, 'AG': 9, 'CA': 9, 'GA': 9, 'TA': 9}
    >>> dict(computingfrequenciesII('TAGCG', 2, 1))
    {'AA': 2, 'AC': 2, 'GT': 1, 'AG': 2, 'CC': 2, 'CA': 2, 'CG': 2, 'GG': 3, 'GC': 1, 'AT': 1, 'GA': 2, 'CT': 1, 'TT': 1, 'TG': 3, 'TC': 2, 'TA': 1}
    >>> dict(computingfrequenciesII('AAAAAAAAAA', 2, 1, True))
    {'AA': 9, 'AC': 9, 'GT': 9, 'AG': 9, 'CA': 9, 'TC': 9, 'AT': 18, 'GA': 9, 'TG': 9, 'TA': 18, 'TT': 9, 'CT': 9}
    """
    frequencyarray = defaultdict(int)
    n = len(text)
    if end is None:
        end = n
    else:
        end = min(n, end)
    for i in range(start, end - k + 1):
        pattern = text[i:i+k]
        frequencyarray[pattern] += 1
    if d > 0: # avoid unneeded step if d == 0
        mfrequencyarray = frequencyarray.copy()
        patterns = frequencyarray.keys()
        for pattern in patterns:
            for neighbor in neighbors(pattern, d):
                if neighbor != pattern:
                    mfrequencyarray[pattern] += frequencyarray[neighbor]
                    if neighbor not in patterns: # if neighbor in patterns, we will count it later, or has been counted before
                        mfrequencyarray[neighbor] += frequencyarray[pattern]
        frequencyarray = mfrequencyarray
    if rc:
        mfrequencyarray = frequencyarray.copy()
        patterns = frequencyarray.keys()
        for pattern in patterns:
            rcpattern = rcomplement.rcomplement(pattern)
            mfrequencyarray[pattern] += frequencyarray[rcpattern]
            if rcpattern not in patterns:
                mfrequencyarray[rcpattern] += frequencyarray[pattern]
        frequencyarray = mfrequencyarray
    return frequencyarray


def frequentwordsII(text, k, d=0, rc=False, start=0, end=None):
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
    >>> frequentwordsII('ACGTTGCATGTCGCATGATGCATGAGAGCT', 4, 1, True)[0]
    set(['ACAT', 'ATGT'])
    >>> frequentwordsII('AAAAAAAAAA', 2, 1, True)[0]
    set(['AT', 'TA'])
    >>> frequentwordsII('AGTCAGTC', 4, 2, True)[0]
    set(['GGCC', 'AATT'])
    >>> frequentwordsII('AATTAATTGGTAGGTAGGTA', 4, 0, True)[0]
    set(['AATT'])
    >>> s = frequentwordsII('ATA', 3, 1, True)[0]
    >>> s == set('AAA AAT ACA AGA ATA ATC ATG ATT CAT CTA GAT GTA TAA TAC TAG TAT TCT TGT TTA TTT'.split(' '))
    True
    >>> frequentwordsII('AAT', 3, 0, True)[0]
    set(['AAT', 'ATT'])
    >>> frequentwordsII('TAGCG', 2, 1, True)[0]
    set(['CC', 'GG', 'CA', 'TG'])
    """
    frequencyarray = computingfrequenciesII(text, k, d, rc, start, end)
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
    parser.add_argument("--hamming", type=int, default=0)
    parser.add_argument("--with-reverse-complement", action='store_true')
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("-L", type=int)
    args = parser.parse_args()

    if args.L is None:
        end = None
    else:
        end = args.start + args.L
    s = multifile.read(args.text)
    fwords = frequentwordsII(s.read().strip(), args.k, args.hamming, args.with_reverse_complement, args.start, end)
    print "frequent words: %s" % ' '.join(fwords[0])
    print "Max frequency: %s" % fwords[1]
