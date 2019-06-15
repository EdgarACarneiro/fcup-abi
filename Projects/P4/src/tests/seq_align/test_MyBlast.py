import unittest

from seq_align import MyBlast


class test_BioSeq(unittest.TestCase):

    def test_create_bio_seq(self):
        pass
        # self.assertIsInstance(BioSeq.create_bio_seq("ATCGTAT"), Dna)
        # self.assertIsInstance(BioSeq.create_bio_seq("ATCGTAT", "dna"), Dna)
        # self.assertIsInstance(BioSeq.create_bio_seq("AUGCUAC", "RNA"), Rna)
        # self.assertIsInstance(BioSeq.create_bio_seq(
        #     "AUNED_MIAEDMA_IAED", "Protein"), Protein)
        # self.assertRaises(BioSeq.BioTypeException, BioSeq.create_bio_seq, "ATCTG", "noseq")

    def test_read_fasta_file(self):
        pass
        # # TP53 AND "Homo sapiens"
        # fasta = BioSeq.read_fasta_file('tests/files/sequence.fasta')
        # self.assertEqual(len(fasta), 15)
        # self.assertEqual(list(fasta.keys())[2], "NM_001126113.2 Homo sapiens tumor protein p53 (TP53), transcript variant 4, mRNA")

    # Test of method load() is made on test_Seq.py

if __name__ == '__main__':
    unittest.main()
