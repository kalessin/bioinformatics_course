"""
Counts the ocurrences of a pattern in a text
"""
import argparse
import hamming

def patterncount(text, pattern, d):
    count = 0
    k = len(pattern)
    n = len(text)
    for i in range(0, n - k + 1):
        if hamming.hamming(text[i:i+k], pattern) <= d:
            count += 1
    return count


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("text", metavar="TEXTFILE", type=argparse.FileType('r'))
    parser.add_argument("pattern")
    parser.add_argument("--hamming", type=int, default=0)

    args = parser.parse_args()
    print "count: %d" % patterncount(args.text.read(), args.pattern, args.hamming)
