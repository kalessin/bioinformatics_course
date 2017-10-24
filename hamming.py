def hamming(seq1, seq2):
    result = 0
    for i, c in enumerate(seq1):
        if c != seq2[i]:
            result += 1
    return result
