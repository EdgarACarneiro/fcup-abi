import unittest

from bioseq import NucleotideChain
from bioseq import BioSeq


class test_NucleotidicChain(unittest.TestCase):

    def test_init(self):
        # Launches Exception since class is abstract
        self.assertRaises(Exception, NucleotideChain, "ACTG")

    def test_validate(self):
        dna = BioSeq.create_bio_seq("ACTG")
        rna = BioSeq.create_bio_seq("ACUG", "rna")
        self.assertTrue(dna.validate())
        self.assertTrue(rna.validate())

        # seq.csv contains a DNA sequence
        rna.read_sequence('tests/files/seq.csv')
        self.assertFalse(rna.validate())

    def test_gc_percent(self):
        dna = BioSeq.create_bio_seq("ACGG")
        self.assertEqual(dna.gc_percent(), 0.75)
        self.assertEqual(dna.gc_percent_sub_seq(2), [0.5, 1])

    def test_rev_complement(self):
        self.assertEqual(BioSeq.create_bio_seq(
            "ACGGTA").reverse_complement(), "TACCGT")
        self.assertEqual(BioSeq.create_bio_seq(
            "ACGUUA", "rna").reverse_complement(), "UAACGU")


if __name__ == '__main__':
    unittest.main()
