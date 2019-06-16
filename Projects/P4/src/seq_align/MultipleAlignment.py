from bioseq import Seq, BioSeq
from seq_align import MyAlign


class MultipleAlignment():
    """Implements the simplified progressive MSA"""

    seqs: [Seq]
    pair_align_data: [list]
    GAP_IDX = 0
    SM_IDX = 1

    def __init__(self, seqs, align_data):
        self.seqs = seqs  # list of MySeq objects
        self.pair_align_data = align_data  # (sm, g) object

    def num_seqs(self):
        """Returns the number of sequnces used in the MSA"""
        return len(self.seqs)

    def add_seq_alignment(self, alignment: MyAlign, seq: Seq):
        """Adds new sequences to existing alignments.
        Aligns the consensus of the previous alignment + new
        sequence = new alignment.
        Returns the new alignment"""
        res = ["" for i in range(len(alignment.get_seqs())+1)]

        # create consensus from given alignments
        cons = BioSeq.create_bio_seq(
            alignment.consensus(), alignment.get_align_type())

        align2 = MyAlign.align_from_global_alignment(
            cons, seq, *self.pair_align_data, alignment.get_align_type())

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
            res = self.add_seq_alignment(res, self.seqs[i])

        return res

    def ScoreColumn(self, charsCol):
        """Calculate the score of each column in the alignment.
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

    def scoreSP(self, alignment: MyAlign):
        """Returns the score SP from a complete alignment"""
        return sum(
            [self.ScoreColumn(alignment.column(col_idx))
             for col_idx in range(len(alignment))]
        )


def test_prot():
    s1 = MySeq("PHWAS", "protein")
    s2 = MySeq("HWASW", "protein")
    s3 = MySeq("HPHWA", "protein")
    sm = SubstMatrix()
    sm.read_submat_file("blosum62.mat", "\t")
    aseq = PairwiseAlignment(sm, -8)
    ma = MultipleAlignment([s1, s2, s3], aseq)
    alinm = ma.align_consensus()
    print(alinm)


def test():
    s1 = MySeq("ATAGC")
    s2 = MySeq("AACC")
    s3 = MySeq("ATGAC")

    sm = SubstMatrix()
    sm.create_submat(1, -1, "ACGT")
    aseq = PairwiseAlignment(sm, -1)
    ma = MultipleAlignment([s1, s2, s3], aseq)
    al = ma.align_consensus()
    print(al)


def exercise1():
    s1 = MySeq("ACATATCAT")
    s2 = MySeq("AACAGATCT")
    s3 = MySeq("AGATATTAG")
    s4 = MySeq("GCATCGATT")

    sm = SubstMatrix()
    sm.create_submat(1, -1, "ACGT")
    aseq = PairwiseAlignment(sm, -1)
    ma = MultipleAlignment([s1, s2, s3, s4], aseq)
    al = ma.align_consensus()
    print(al)


if __name__ == "__main__":
    test_prot()
    print()
    test()
    print()
    exercise1()
