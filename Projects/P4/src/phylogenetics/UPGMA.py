from phylogenetics import NumMatrix, HierarchicalClustering
from bioseq import Seq
from seq_align import SubstMatrix, MyAlign


class UPGMA:

    seqs: [Seq]
    align_data: (SubstMatrix, float)
    dists_mat: NumMatrix

    def __init__(self, seqs, align_data):
        self.seqs = seqs
        self.align_data = align_data
        self.__create_mat_dist()

    def __create_mat_dist(self):
        # create distance matrix
        self.dists_mat = NumMatrix(self.seqs, len(self.seqs))

        for i in range(len(self.seqs)):
            for j in range(i, len(self.seqs)):
                # align the two sequences and recover the alignment
                align = MyAlign.align_from_global_alignment(
                    self.seqs[i], self.seqs[j], *self.align_data)

                num_different_symbols = 0

                # count the number of different symbols in the alignment as the distance
                for k in range(len(align)):
                    col = align.column(k)
                    if (col[0] != col[1]):
                        num_different_symbols += 1

                # set distance value in the matrix
                self.dists_mat.set_value(i, j, num_different_symbols)

    def run(self):
        # create an object of the class HierarchicalClustering
        # with the computed distance matrix
        cluster = HierarchicalClustering(self.dists_mat)

        # execute the clustering algorithm
        return cluster.execute_clustering()


def test():
    seq1 = MySeq("ATAGCGAT")
    seq2 = MySeq("ATAGGCCT")
    seq3 = MySeq("CTAGGCCC")
    seq4 = MySeq("CTAGGCCT")
    sm = SubstMatrix()
    sm.create_submat(1, -1, "ACGT")
    align_data = PairwiseAlignment(sm, -2)
    up = UPGMA([seq1, seq2, seq3, seq4], align_data)
    arv = up.run()
    arv.print_tree()


# exercise 1c) of chapter 9
def exercise1():
    s1 = MySeq("ACATATCAT")
    s2 = MySeq("AACAGATCT")
    s3 = MySeq("AGATATTAG")
    s4 = MySeq("GCATCGATT")

    sm = SubstMatrix()
    sm.create_submat(1, -1, "ACGT")
    aseq = PairwiseAlignment(sm, -1)

    up = UPGMA([s1, s2, s3, s4], aseq)
    arv = up.run()
    arv.print_tree()


if __name__ == '__main__':
    test()
    print()
    # exercise1()
