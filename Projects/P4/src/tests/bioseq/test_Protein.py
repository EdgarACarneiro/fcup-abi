import unittest

from bioseq import Protein
from bioseq import Seq


class test_Protein(unittest.TestCase):

    def test_init(self):
        self.assertEqual(Protein("_INAED").get_seq(), "Protein: _INAED")
        # Launches Exception since sequence contains invalid chars ('3'/ '!')
        self.assertRaises(Seq.InvalidSequenceException, Protein, "ACTG32")
        self.assertRaises(Seq.InvalidSequenceException, Protein, "ACSG!")

    def test_basics(self):
        self.assertEqual(str(Protein("_INAED")), "_INAED")

    def test_rev_transcription(self):
        self.assertEqual(Protein("INAEDMADME_").all_proteins_rf(), ["MADME", "ME"])

if __name__ == '__main__':
    unittest.main()
