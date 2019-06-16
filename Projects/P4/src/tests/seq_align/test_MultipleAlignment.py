import unittest

from bioseq import Dna, Rna, Protein
from seq_align import MyAlign, SubstMatrix, MultipleAlignment


class test_MultipleAlignment(unittest.TestCase):

    def test_align_consensus(self):
        """Test that iteratively tests other functions"""
        s1 = Protein("PHWAS")
        s2 = Protein("HWASW")
        s3 = Protein("HPHWA")
        sm = SubstMatrix.read_submat_file("tests/files/blosum62.mat")
        ma = MultipleAlignment([s1, s2, s3], (sm, -8))
        self.assertEqual(
            str(ma.align_consensus()),
            '\n-PHWAS-\n--HWASW\nHPHWA--'
        )

        s1 = Dna("ATAGC")
        s2 = Dna("AACC")
        s3 = Dna("ATGAC")
        sm = SubstMatrix.create_submat("ACGT", 1, -1)
        ma = MultipleAlignment([s1, s2, s3], (sm, -1))
        self.assertEqual(
            str(ma.align_consensus()),
            '\nAT-AGC\nA--ACC\nATGA-C'
        )

        s1 = Dna("ACATATCAT")
        s2 = Dna("AACAGATCT")
        s3 = Dna("AGATATTAG")
        s4 = Dna("GCATCGATT")
        ma = MultipleAlignment([s1, s2, s3, s4], (sm, -1))
        self.assertEqual(
            str(ma.align_consensus()),
            '\n-AC-AT--ATCAT\nAAC-AG--ATC-T\n-AG-AT--ATTAG\n--GCATCGA-T-T'
        )

        print('>> Passed test_align_consensus()')

    def test_score_sp(self):
        alignment = MyAlign([Protein("PHWAS"), Protein("HWASW"), Protein("HPHWA")])
        sm = SubstMatrix.read_submat_file("tests/files/blosum62.mat")
        ma = MultipleAlignment([], (sm, -8))

        self.assertEqual(ma.score_sp(alignment), -21)

        print('>> Passed test_score_sp()')


if __name__ == '__main__':
    unittest.main()
