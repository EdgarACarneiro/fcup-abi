from functools import reduce

def read_submat_file(filename):
    """Read a substitution matrix from the given file"""
    f = open(filename, "r")
    alphabet = [symbol for symbol in f.readline().replace('\n', '').split('\t')]

    dic = {}
    for i, line in enumerate(f):
        line_symbol = line.replace('\n', '').split('\t')
        
        for j in range(0, len(line_symbol)):
            dic[alphabet[i] + alphabet[j]] = int(line_symbol[j])

    f.close()
    return dic


def __max3(v1, v2, v3):
    """Indicates, between the given integers, the indexes of the biggest ones"""
    if v1 > v2:
        return 1 if v1 > v3 else 3
    else:
        return 2 if v2 > v3 else 3


def __score_pos(c1, c2, sm, g):
    """Score of a column alignment (between c1 and c2).
    Assume a constant gap penalty g and a substituin matrix sm"""
    return g if c1 == '-' or c2 == '-' else sm[c1+c2]


def __score_align(seq1, seq2, sm, g):
    """Score of the whole sequence alignment"""
    assert len(seq1) == len(seq2), "Sequences must have same size"
    
    return reduce((lambda acc, val: acc + val), \
                  [__score_pos(seq1[i], seq2[i], sm, g) for i in range(0, len(seq1))])


def global_align_multiple_solutions(seq1, seq2, sm, g):
    """Global Alignment with multiple solutions of two sequences"""
    score = [[0]]
    trace = [[0]]

    # initialize gaps in cols
    for i in range(1, len(seq1) + 1):
        score.append([g * i])
        trace.append([2])
    
    # initialize gaps in rows
    for j in range(1, len(seq2) + 1):
        score[0].append(g * j)
        trace[0].append(3)

    # apply the recurrence to fill the matrices
    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            v1 = score[i-1][j-1] + __score_align(seq1[i-1], seq2[j-1], sm, g)
            v2 = score[i-1][j] + g
            v3 = score[i][j-1] + g
            score[i].append(max(v1, v2, v3))
            trace[i].append(__max3(v1, v2, v3))

    return (score, trace)