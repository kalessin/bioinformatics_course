import argparse

import multifile


def skew(text):
    skew = 0
    yield skew
    for c in text:
        if c == 'C':
            skew -= 1
        elif c == 'G':
            skew += 1
        yield skew


def minskewindex(text, reverse=False):
    result = []
    minskew = float('inf')
    n = len(text)
    if reverse:
        text = text[::-1]
    for i, s in enumerate(skew(text)):
        if s < minskew:
            result = [i]
            minskew = s
        elif s == minskew:
            result.append(i)
    if reverse:
        result = [n - i for i in result]
    return result


def maxskewindex(text, reverse=False):
    result = []
    maxskew = -float('inf')
    n = len(text)
    if reverse:
        text = text[::-1]
    for i, s in enumerate(skew(text)):
        if s > maxskew:
            result = [i]
            maxskew = s
        elif s == maxskew:
            result.append(i)
    if reverse:
        result = [n - i for i in result]
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('text', type=argparse.FileType('r'))
    parser.add_argument('--reverse', action='store_true')
    parser.add_argument('--find-max', action='store_true')
    args = parser.parse_args()
    s = multifile.read(args.text)
    if args.find_max:
        print ' '.join(`i` for i in maxskewindex(s.read().strip(), args.reverse))
    else:
        print ' '.join(`i` for i in minskewindex(s.read().strip(), args.reverse))
