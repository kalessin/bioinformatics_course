import argparse

def find(text, pattern):
    pos = []
    n = len(text)
    k = len(pattern)
    for i in range(0, n - k + 1):
        if text[i:i+k] == pattern:
            pos.append(i)
    return pos


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("text", metavar="TEXTFILE", type=argparse.FileType('r'))
    parser.add_argument("pattern")

    args = parser.parse_args()
    print ' '.join(`i` for i in find(args.text.read().strip(), args.pattern))
