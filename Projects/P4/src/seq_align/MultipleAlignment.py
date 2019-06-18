from ..bioseq import Seq, BioSeq
from ..seq_align import MyAlign, SubstMatrix


class MultipleAlignment():
    """Implements the simplified progressive MSA"""

    seqs: [Seq]
    pair_align_data: (SubstMatrix, float)
    GAP_IDX = 1
    SM_IDX = 0

    def __init__(self, seqs, align_data):
        self.seqs = seqs  # list of MySeq objects
        self.pair_align_data = align_data  # (sm, g) object

    def num_seqs(self):
        """Returns the number of sequnces used in the MSA"""
        return len(self.seqs)

    def __add_seq_alignment(self, alignment: MyAlign, seq: Seq):
        """Adds new sequences to existing alignments.
        Aligns the consensus of the previous alignment + new
        sequence = new alignment.
        Returns the new alignment"""
        res = ["" for i in range(len(alignment.get_seqs())+1)]

        # create consensus from given alignments
        cons = BioSeq.create_bio_seq(
            alignment.consensus(), alignment.get_align_type())

        align2 = MyAlign.align_from_global_alignment(
            cons, seq, *self.pair_align_data)

        orig = 0
        for i in range(len(align2)):
            if align2[0, i] == "-":
                for k in range(len(alignment.get_seqs())):
                    res[k] += "-"
            else:
                for k in range(len(alignment.get_seqs())):
                    res[k] += alignment[k, orig]
                orig += 1

        res[len(alignment.get_seqs())] = align2.get_seqs()[1]

        return MyAlign(res, alignment.get_align_type())

    def align_consensus(self):
        """General implementation of the MSA.
        Progressive consensus on the sequence list."""
        res = MyAlign.align_from_global_alignment(
            self.seqs[0], self.seqs[1], *self.pair_align_data)

        for i in range(2, len(self.seqs)):
            res = self.__add_seq_alignment(res, self.seqs[i])

        return res

    def __score_column(self, charsCol):
        """Calculate the score of the given column.
        CharsCol - list of chars representing a column in an alignment.
        The score is computed using the Sum of Pairs (SP) approach,
        i.e.the score will be the sum of the scores of each pair of
        characters in the alignment. If two gaps are found in each
        pair then the score will be zero."""
        score = 0

        for i in range(0, len(charsCol)-1):
            for j in range(i+1, len(charsCol)):
                p1 = charsCol[i]
                p2 = charsCol[j]

                if p1 == '-' and p2 == '-':
                    pass
                elif p1 == '-' or p2 == '-':
                    score += self.pair_align_data[self.GAP_IDX]
                else:
                    score += self.pair_align_data[self.SM_IDX][p1, p2]

        return score

    def score_sp(self, alignment):
        """Returns the score SP from each column of a complete alignment"""
        return sum(
            [self.__score_column(alignment.column(col_idx))
             for col_idx in range(len(alignment))]
        )
