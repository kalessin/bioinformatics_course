import argparse
from collections import defaultdict


def computingfrequenciesII(text, k):
    """
    Similar to frequentwords.computingfrequenciesII, but instead of counting, it generates
    a list containing the occurrences indexes
    """
    frequencyarray = defaultdict(list)
    n = len(text)
    for i in range(0, n - k + 1):
        pattern = text[i:i+k]
        frequencyarray[pattern].append(i)
    return frequencyarray


def clumpfind(text, k, t, L):
    frequencyarray = computingfrequenciesII(text, k)
    clumps = defaultdict(list)
    for _p, _o in frequencyarray.iteritems():
        if len(_o) < t:
            continue
        for i in range(0, len(_o) - t + 1):
            if _o[i+t-1] + k - _o[i] <= L:
                clumps[_p].append(_o[i:i+t])

    return clumps


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("text", metavar="TEXTFILE", type=argparse.FileType('r'))
    parser.add_argument("k", type=int)
    parser.add_argument("t", type=int)
    parser.add_argument("L", type=int)

    args = parser.parse_args()
    result = clumpfind(args.text.read().strip(), args.k, args.t, args.L)
    candidates = result.keys()
    print "Number of candidates:", len(candidates)
    print ' '.join(candidates)
