import unittest

from ...bioseq import Dna
from ...seq_align import SubstMatrix
from ...phylogenetics import UPGMA


class test_UPGMA(unittest.TestCase):

    def test_init_UPGMA(self):
        seq1 = Dna("ATAGCGAT")
        seq2 = Dna("ATAGGCCT")
        seq3 = Dna("CTAGGCCC")
        seq4 = Dna("CTAGGCCT")
        sm = SubstMatrix.create_submat("ACGT", 1, -1)

        up = UPGMA([seq1, seq2, seq3, seq4], (sm, -2))
        self.assertTrue(seq4 in up.seqs)
        self.assertEqual(up.align_data, (sm, -2))
        self.assertEqual(up.dists_mat[3], [4, 1, 1, 0])

        print('>> Passed test_init_UPGMA()')

    def test_run_UPGMA(self):
        seq1 = Dna("ATAGCGAT")
        seq2 = Dna("ATAGGCCT")
        seq3 = Dna("CTAGGCCC")
        seq4 = Dna("CTAGGCCT")
        sm = SubstMatrix.create_submat("ACGT", 1, -1)
        tree = UPGMA([seq1, seq2, seq3, seq4], (sm, -2)).run()
        self.assertTrue(tree.exists_leaf(1))
        self.assertTrue(tree.exists_leaf(2))
        self.assertTrue(tree.exists_leaf(3))
        self.assertEqual(tree.distance_leaves(3, 0), 4.0)

        s1 = Dna("ACATATCAT")
        s2 = Dna("AACAGATCT")
        s3 = Dna("AGATATTAG")
        s4 = Dna("GCATCGATT")
        sm = SubstMatrix.create_submat("ACGT", 1, -1)
        tree = UPGMA([s1, s2, s3, s4], (sm, -1)).run()
        self.assertTrue(tree.exists_leaf(1))
        self.assertTrue(tree.exists_leaf(2))
        self.assertTrue(tree.exists_leaf(3))
        self.assertEqual(tree.distance_leaves(2, 0), 4.5)

        print('>> Passed test_run_UPGMA()')


if __name__ == '__main__':
    unittest.main()
