from functools import reduce


def __matrix(row, col, val=0):
    """ create a matrix with the given dimensions: 
    number of rows and columns filled with the given value (by omission zero) """
    return [[val] * col for _ in range(0, row)]


def read_submat_file(filename):
    """Read a substitution matrix from the given file"""
    f = open(filename, "r")
    alphabet = [symbol for symbol in f.readline().replace('\n',
                                                          '').split('\t')]

    dic = {}
    for i, line in enumerate(f):
        line_symbol = line.replace('\n', '').split('\t')

        for j in range(0, len(line_symbol)):
            dic[alphabet[i] + alphabet[j]] = int(line_symbol[j])

    f.close()
    return dic


def subst_matrix(alphabet, match, mismatch):
    """Substitution matrix as a dictionary"""
    return {i + j: match if i == j else mismatch
            for i in alphabet for j in alphabet}


def pretty_matrix(matrix, row_label, col_label):
    """Pretty print of the given matrix """

    # Restraining labels that are too big
    row_label = [el[:10] + '..' if len(el) > 10 else el
                 for el in row_label]
    col_label = [el[:10] + '..' if len(el) > 10 else el
                 for el in col_label]

    # Stringfying everything & Joining top label
    s_matrix = [list([" "] + (col_label))] + \
               [[row_label[row_idx]] +
                [str(e) for e in row] for row_idx, row in enumerate(matrix)]

    # Length of each matrix column
    len_s = [max(map(len, col)) for col in zip(*s_matrix)]

    # Cell formatation
    formatation = '\t'.join('{{:{}}}'.format(x) for x in len_s)

    # Apply cell formation to each matrix element
    pretty_mat = [formatation.format(*row) for row in s_matrix]

    # Print Pretty Matrix
    print('\n'.join(pretty_mat))


def __max3(v1, v2, v3):
    """Indicates, between the given integers, the indexes of the biggest ones"""
    vals = [v1, v2, v3]

    return [index + 1 for index, value in enumerate(vals) if value == max(vals)]


def __score_pos(c1, c2, sm, g):
    """Score of a column alignment (between c1 and c2).
    Assume a constant gap penalty g and a substituin matrix sm"""
    return g if c1 == '-' or c2 == '-' else sm[c1+c2]


def __score_align(seq1, seq2, sm, g):
    """Score of the whole sequence alignment"""
    assert len(seq1) == len(seq2), "Sequences must have same size"

    return reduce((lambda acc, val: acc + val),
                  [__score_pos(seq1[i], seq2[i], sm, g) for i in range(0, len(seq1))])


def global_align_multiple_solutions(seq1, seq2, sm, g):
    """Global Alignment with multiple solutions of two sequences.
    Represents an adaptation of the needleman-wunsch algorithm"""
    score = [[0]]
    trace = [[[0]]]  # Each element of the trace matrix is now a list

    # initialize gaps in cols
    for i in range(1, len(seq1) + 1):
        score.append([g * i])
        trace.append([[2]])

    # initialize gaps in rows
    for j in range(1, len(seq2) + 1):
        score[0].append(g * j)
        trace[0].append([3])

    # apply the recurrence to fill the matrices
    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            v1 = score[i-1][j-1] + __score_align(seq1[i-1], seq2[j-1], sm, g)
            v2 = score[i-1][j] + g
            v3 = score[i][j-1] + g
            score[i].append(max(v1, v2, v3))
            trace[i].append(__max3(v1, v2, v3))

    return (score, trace)


__DIAGONAL = 1
__VERTICAL = 2
__HORIZONTAL = 3


def __recover_global_aux(trace, seq1, seq2, memoization_mat):
    """Auxiliary recursive function with MEMOIZATION to help return all the 
    possible global align solutions"""
    i, j = len(seq1), len(seq2)
    res = []

    if (i > 0 or j > 0):
        if len(memoization_mat[i][j]) == 0:
            if __DIAGONAL in trace[i][j]:
                res += [
                    [val[0] + seq1[i-1], val[1] + seq2[j-1]]
                    for val
                    in __recover_global_aux(trace, seq1[:-1], seq2[:-1], memoization_mat)
                ]

            if __VERTICAL in trace[i][j]:
                res += [
                    [val[0] + seq1[i-1], val[1] + '-']
                    for val
                    in __recover_global_aux(trace, seq1[:-1], seq2, memoization_mat)
                ]

            if __HORIZONTAL in trace[i][j]:
                res += [
                    [val[0] + '-', val[1] + seq2[j-1]]
                    for val
                    in __recover_global_aux(trace, seq1, seq2[:-1], memoization_mat)
                ]

            # Updating the memoization matrix with the computed arrays
            memoization_mat[i][j] = res

        return memoization_mat[i][j]

    # Base case position 0, 0
    return [["", ""]]


def recover_global_align_multiple_solutions(trace, seq1, seq2):
    """Recover the optimal alignments between seq1 and seq2
    using the given trace matrix and recursivity"""
    return __recover_global_aux(trace, seq1, seq2,
                                __matrix(len(seq1) + 1, len(seq2) + 2, []))


def local_align_multiple_solutions(seq1, seq2, sm, g):
    """Local Alignment with multiple solutions of two sequences.
    Represents an adaptation of the smith-waterman algorithm"""
    score = [[0]]
    trace = [[[0]]]  # Each element of the trace matrix is now a list

    # initialize gaps in cols
    for i in range(1, len(seq1) + 1):
        score.append([0])
        trace.append([[0]])
        max_score = 0

    # initialize gaps in rows
    for j in range(1, len(seq2) + 1):
        score[0].append(0)
        trace[0].append([0])

    # apply the recurrence to fill the matrices
    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            v1 = score[i-1][j-1] + __score_align(seq1[i-1], seq2[j-1], sm, g)
            v2 = score[i-1][j] + g
            v3 = score[i][j-1] + g
            max_v = max(v1, v2, v3)

            if (max_v <= 0):
                score[i].append(0)
                trace[i].append([0])
            else:
                score[i].append(max_v)
                trace[i].append(__max3(v1, v2, v3))

            max_score = max_score if max_v < max_score else max_v

    return (score, trace, max_score)


def __max_positions(mat):
    """Returns all the positions of the existing maximum value in
    the given matrix"""
    pos = []
    max_val = mat[0][0]

    for i in range(0, len(mat)):
        for j in range(0, len(mat[i])):
            if mat[i][j] > max_val:
                max_val = mat[i][j]
                pos = [(i, j)]

            elif mat[i][j] == max_val:
                pos.append((i, j))

    return (pos, max_val)


def __recover_local_aux(trace, seq1, seq2, pos, memoization_mat):
    """Auxiliary recursive function with MEMOIZATION to help return all the 
    possible global align solutions"""
    i, j = pos
    res = []

    if (trace[i][j] > [0]):
        if len(memoization_mat[i][j]) == 0:
            if __DIAGONAL in trace[i][j]:
                res += [
                    [val[0] + seq1[i-1], val[1] + seq2[j-1]]
                    for val
                    in __recover_local_aux(trace, seq1, seq2, (i-1, j-1), memoization_mat)
                ]

            if __VERTICAL in trace[i][j]:
                res += [
                    [val[0] + seq1[i-1], val[1] + '-']
                    for val
                    in __recover_local_aux(trace, seq1, seq2, (i-1, j), memoization_mat)
                ]

            if __HORIZONTAL in trace[i][j]:
                res += [
                    [val[0] + '-', val[1] + seq2[j-1]]
                    for val
                    in __recover_local_aux(trace, seq1, seq2, (i, j-1), memoization_mat)
                ]

            # Updating the memoization matrix with the computed arrays
            memoization_mat[i][j] = res

        return memoization_mat[i][j]

    # Base case position trace[i][j] == 0
    return [["", ""]]


def recover_local_align_multiple_solutions(score, trace, seq1, seq2):
    la_pos, _ = __max_positions(score)  # Local Alignment starting positions
    memoization_mat = __matrix(len(seq1) + 1, len(seq2) + 2, [])

    return reduce((lambda acc, val: acc + val),
                  [__recover_local_aux(trace, seq1, seq2, pos, memoization_mat)
                   for pos in la_pos])


def compare_pairwise_global_align(seq_list, sm, g):
    """Gets a matrix indicating the maximum score of possible global alignments
    resultant of the cross product of sequences"""
    cross_prod = []

    for i in range(0, len(seq_list)):
        cross_prod.append([])

        for seq2 in seq_list:
            ga_score, _ = global_align_multiple_solutions(
                seq_list[i], seq2, sm, g)
            _, max_val = __max_positions(ga_score)

            cross_prod[i].append(max_val)

    pretty_matrix(cross_prod, seq_list, seq_list)

    return cross_prod


def compare_pairwise_local_align(seq_list, sm, g):
    """Gets a matrix indicating the score of possible local alignments
    resultant of the cross product of sequences"""
    cross_prod = []

    for i in range(0, len(seq_list)):
        cross_prod.append([])

        for seq2 in seq_list:
            *_, max_val = local_align_multiple_solutions(seq_list[i], seq2, sm, g)
            cross_prod[i].append(max_val)

    pretty_matrix(cross_prod, seq_list, seq_list)

    return cross_prod


def compare_pairwise_num_global_align(seq_list, sm, g):
    """Gets a matrix indicating the number of possible global alignments
    resultant of the cross product of sequences"""
    cross_prod = []

    for i in range(0, len(seq_list)):
        cross_prod.append([])

        for seq2 in seq_list:
            _, ga_trace = global_align_multiple_solutions(
                seq_list[i], seq2, sm, g)
            cross_prod[i].append(
                len(recover_global_align_multiple_solutions(
                    ga_trace, seq_list[i], seq2))
            )

    pretty_matrix(cross_prod, seq_list, seq_list)

    return cross_prod


def compare_pairwise_num_local_align(seq_list, sm, g):
    """Gets a matrix indicating the number of possible local alignments
    resultant of the cross product of sequences"""
    cross_prod = []

    for i in range(0, len(seq_list)):
        cross_prod.append([])

        for seq2 in seq_list:
            ga_score, ga_trace, _ = local_align_multiple_solutions(
                seq_list[i], seq2, sm, g)
            cross_prod[i].append(
                len(recover_local_align_multiple_solutions(
                    ga_score, ga_trace, seq_list[i], seq2))
            )

    pretty_matrix(cross_prod, seq_list, seq_list)

    return cross_prod


def align_query(query_seq, seq_list, sm, g):
    """Finds the most similar sequence to a query sequence.
    Uses local alignment"""
    max_result = (-1, None, None)  # score, seq, alignment

    for seq in seq_list:
        la = local_align_multiple_solutions(query_seq, seq, sm, g)
        if la[2] > max_result[0]:
            max_result = (la[2], seq, la)

    best_alignment = recover_local_align_multiple_solutions(
        max_result[2][0], max_result[2][1], query_seq, max_result[1])

    # Return the first of the best alignments
    return best_alignment[0], max_result[0]
