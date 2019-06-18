import unittest

from ...bioseq import Rna
from ...bioseq import Seq
from ...bioseq import Dna


class test_Rna(unittest.TestCase):

    def test_init(self):
        self.assertEqual(Rna("UCGA").get_seq(), "RNA: UCGA")
        # Launches Exception since sequence contains invalid chars ('T'/ 'S')
        self.assertRaises(Seq.InvalidSequenceException, Rna, "ACTG")
        self.assertRaises(Seq.InvalidSequenceException, Rna, "ACSG")

        print('>> Passed Rna::test_init()')

    def test_basics(self):
        self.assertEqual(str(Rna("UCAUU")), "UCAUU")

        print('>> Passed Rna::test_basics()')

    def test_rev_transcription(self):
        dna = Rna("UCAUU").rev_transcription()
        self.assertIsInstance(dna, Dna)
        self.assertEqual(dna.get_seq(), "DNA: TCATT")

        print('>> Passed Rna::test_rev_transcription()')

if __name__ == '__main__':
    unittest.main()
