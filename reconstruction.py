def composition(k, text):
    """
    >>> composition(5, 'CAATCCAAC')
    ['AATCC', 'ATCCA', 'CAATC', 'CCAAC', 'TCCAA']
    """
    n = len(text)
    result = set()
    for i in xrange(n - k + 1):
        result.add(text[i:i+k])
    return sorted(result)


def reconstruct_from_sorted_kmers(kmers):
    """
    >>> reconstruct_from_sorted_kmers(['ACCGA', 'CCGAA', 'CGAAG', 'GAAGC', 'AAGCT'])
    'ACCGAAGCT'
    """
    result, kmers = kmers[0], kmers[1:]
    for kmer in kmers:
        k = len(kmer)
        assert kmer[0:k-1] == result[-k+1:]
        result += kmer[-1]
    return result
