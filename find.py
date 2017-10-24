"""
Finds the position of a pattern in a text
"""
import argparse

import hamming

def find(text, pattern, hamming_distance=0):
    pos = []
    n = len(text)
    k = len(pattern)
    for i in range(0, n - k + 1):
        if hamming.hamming(text[i:i+k], pattern) <= hamming_distance:
            pos.append(i)
    return pos


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("text", metavar="TEXTFILE", type=argparse.FileType('r'))
    parser.add_argument("pattern")
    parser.add_argument("--hamming", type=int, default=0)

    args = parser.parse_args()
    print ' '.join(`i` for i in find(args.text.read().strip(), args.pattern, args.hamming))
