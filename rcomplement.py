import argparse

_MAP = ('ATGC', 'TACG')

def scomplement(s):
    return _MAP[1][_MAP[0].find(s)]

def rcomplement(pattern):
    return ''.join(scomplement(s) for s in pattern[::-1])


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--text", metavar="TEXTFILE", type=argparse.FileType('r'))
    parser.add_argument("--pattern")

    args = parser.parse_args()
    print rcomplement(args.text.read().strip() if args.text else args.pattern)
