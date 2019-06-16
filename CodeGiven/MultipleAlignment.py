from PairwiseAlignment import PairwiseAlignment
from MyAlign import MyAlign
from MySeq import MySeq
from SubstMatrix import SubstMatrix


class MultipleAlignment():

    def __init__(self, seqs, alignseq):
        self.seqs = seqs  # list of MySeq objects
        self.alignpars = alignseq  # PairwiseAlignment objects

    def num_seqs(self):
        return len(self.seqs)

    def add_seq_alignment(self, alignment: MyAlign, seq):
        """Adds new sequences to existing alignments. Used as: Aligns
        the consensus of the previous alignment + new sequence = new alignemnt"""
        res = []
        for i in range(len(alignment.listseqs)+1):
            res.append("")
        # create consensus from given alignments
        cons = MySeq(alignment.consensus(), alignment.al_type)
        self.alignpars.needleman_Wunsch(cons, seq)
        align2 = self.alignpars.recover_align()
        orig = 0
        for i in range(len(align2)):
            if align2[0, i] == "-":
                for k in range(len(alignment.listseqs)):
                    res[k] += "-"
            else:
                for k in range(len(alignment.listseqs)):
                    res[k] += alignment[k, orig]
                orig += 1
        res[len(alignment.listseqs)] = align2.listseqs[1]
        return MyAlign(res, alignment.al_type)

    def align_consensus(self):
        """General implmentation of the MSA."""
        self.alignpars.needleman_Wunsch(self.seqs[0], self.seqs[1])
        res = self.alignpars.recover_align()

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
                    score += self.alignpars.g
                else:
                    score += self.alignpars.sm[p1, p2]
        return score

    def scoreSP(self, alignment):
        """Returns the score SP from a complete alignment"""
        # for el in alignment:
        #     self.SCoreColumn()
        # TODO
        return None


def printMat(mat):
    for i in range(0, len(mat)):
        print(mat[i])


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
