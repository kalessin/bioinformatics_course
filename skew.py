import argparse


def skew(text):
    skew = 0
    yield skew
    for c in text:
        if c == 'C':
            skew -= 1
        elif c == 'G':
            skew += 1
        yield skew

def minskewindex(text):
    result = []
    minskew = float('inf')
    for i, s in enumerate(skew(text)):
        if s < minskew:
            result = [i]
            minskew = s
        elif s == minskew:
            result.append(i)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('text', type=argparse.FileType('r'))
    args = parser.parse_args()
    print ' '.join(`i` for i in minskewindex(args.text.read().strip()))
