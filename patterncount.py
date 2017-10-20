import argparse

def patterncount(text, pattern):
    count = 0
    k = len(pattern)
    n = len(text)
    for i in range(0, n - k + 1):
        if text[i:i+k] == pattern:
            count += 1
    return count


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("text", metavar="TEXTFILE", type=argparse.FileType('r'))
    parser.add_argument("pattern")

    args = parser.parse_args()
    print "count: %d" % patterncount(args.text.read(), args.pattern)
