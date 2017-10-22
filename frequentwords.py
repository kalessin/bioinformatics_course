"""
Returns the most frequent k-mer patterns in a text, along with the count of ocurrencies
"""
import argparse
from itertools import product
from collections import defaultdict


import patterncount


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


def computingfrequenciesII(text, k):
    frequencyarray = defaultdict(int)
    n = len(text)
    for i in range(0, n - k + 1):
        pattern = text[i:i+k]
        frequencyarray[pattern] += 1
    return frequencyarray


def frequentwordsII(text, k):
    frequencyarray = computingfrequenciesII(text, k)
    maxcount = max(frequencyarray.itervalues())
    if maxcount == 1:
        return set(), 1
    frequentpatterns = set()
    for p, c in frequencyarray.iteritems():
        if c == maxcount:
            frequentpatterns.add(p)
    return frequentpatterns, maxcount

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("text", metavar="TEXTFILE", type=argparse.FileType("r"))
    parser.add_argument("k", type=int)
    parser.add_argument("--alg", type=int, default=0)
    args = parser.parse_args()
    algmap = {
        0: frequentwords,
        1: fasterfrequentwords,
        2: frequentwordsII,
    }

    print "frequent words: %s (%d)" % algmap[args.alg](args.text.read().strip(), args.k)
