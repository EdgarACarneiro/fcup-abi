from phylogenetics import NumMatrix, HierarchicalClustering
from bioseq import Seq
from seq_align import SubstMatrix, MyAlign


class UPGMA:
    """Implementation of the 'Unweighted Pair Group Method
    Using Arithmetic Averages' method"""

    seqs: [Seq]
    align_data: (SubstMatrix, float)
    dists_mat: NumMatrix

    def __init__(self, seqs, align_data):
        self.seqs = seqs
        self.align_data = align_data
        self.__create_mat_dist()

    def __create_mat_dist(self):
        # create distance matrix
        self.dists_mat = NumMatrix(len(self.seqs), len(self.seqs))

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
