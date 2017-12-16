def decomposition(k, text):
    """
    decompose a string into a sequence of overlapping k-mers
    >>> decomposition(5, 'CAATCCAAC')
    ['AATCC', 'ATCCA', 'CAATC', 'CCAAC', 'TCCAA']
    """
    n = len(text)
    result = set()
    for i in xrange(n - k + 1):
        result.add(text[i:i+k])
    return sorted(result)


def compose_from_sorted_kmers(kmers):
    """
    compose original string from overlapping k-mers decomposition sequence
    >>> compose_from_sorted_kmers(['ACCGA', 'CCGAA', 'CGAAG', 'GAAGC', 'AAGCT'])
    'ACCGAAGCT'
    """
    result, kmers = kmers[0], kmers[1:]
    for kmer in kmers:
        k = len(kmer)
        assert kmer[0:k-1] == result[-k+1:]
        result += kmer[-1]
    return result


def overlap_graph(kmers):
    """
    For each k-mer in a sequence, generate a list of overlapping k-mers from it
    >>> overlap_graph(['ATGCG', 'GCATG', 'CATGC', 'AGGCA', 'GGCAT'])
    [('GCATG', ['CATGC']), ('CATGC', ['ATGCG']), ('AGGCA', ['GGCAT']), ('GGCAT', ['GCATG'])]
    >>> overlap_graph(['AAAAA', 'AAAAA', 'AAAAC', 'AAACG', 'AACGA', 'ACGAT'])
    [('AAAAA', ['AAAAA', 'AAAAC']), ('AAAAA', ['AAAAA', 'AAAAC']), ('AAAAC', ['AAACG']), ('AAACG', ['AACGA']), ('AACGA', ['ACGAT'])]
    """
    result = []
    for i in range(len(kmers)):
        kmer = kmers[i]
        k = len(kmer)
        ovp = []
        for j in range(len(kmers)):
            if i != j:
                kkmer = kmers[j]
                if kkmer[0:k-1] == kmer[-k+1:]:
                    ovp.append(kkmer)
        if ovp:
            result.append((kmer, ovp))
    return result
