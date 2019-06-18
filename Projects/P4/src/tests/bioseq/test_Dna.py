import unittest

from ...bioseq import Rna
from ...bioseq import Seq
from ...bioseq import Dna


class test_Rna(unittest.TestCase):

    def test_init(self):
        self.assertEqual(Dna("ATCGAT").get_seq(), "DNA: ATCGAT")
        # Launches Exception since sequence contains invalid chars ('U'/ 'S')
        self.assertRaises(Seq.InvalidSequenceException, Dna, "ATCUGT")
        self.assertRaises(Seq.InvalidSequenceException, Dna, "ACSG")

        print('>> Passed Dna::test_init()')

    def test_basics(self):
        self.assertEqual(str(Dna("ATCGAT")), "ATCGAT")

        print('>> Passed Dna::test_basics()')

    def test_transcription(self):
        rna = Dna("TCATT").transcription()
        self.assertIsInstance(rna, Rna)
        self.assertEqual(rna.get_seq(), "RNA: UCAUU")

        print('>> Passed test_transcription()')

    def test_genetic_code(self):
        dna = Dna("ACTG")
        dna.read_genetic_code('src/tests/files/genetic_code.txt')
        self.assertEqual(len(dna.get_genetic_code()), 64)
        self.assertEqual(list(dna.get_genetic_code())[5], "TGC")
        dna.set_genetic_code({"ABC": "C"})
        self.assertEqual(list(dna.get_genetic_code())[0], "ABC")

        print('>> Passed test_genetic_code()')

    def test_translate(self):
        dna = Dna("ACTGAATCG")
        dna.read_genetic_code('src/tests/files/genetic_code.txt')
        self.assertEqual(dna.translate().get_seq(), "Protein: TES")

        print('>> Passed test_translate()')

    def test_codon_usage(self):
        dna = Dna("AGCTAGCTAGCTACATCAAACGATCGTCGCTAGCTAGCTAC")
        dna.read_genetic_code('src/tests/files/genetic_code.txt')
        self.assertEqual(dna.codon_usage("R"), {'CGT': 0.5, 'CGC': 0.5})

        print('>> Passed test_codon_usage()')

    def test_reading_frames(self):
        dna = Dna("ACGTACGATATGTA")
        dna.read_genetic_code('src/tests/files/genetic_code.txt')
        rf = dna.reading_frames()
        self.assertEqual(len(rf), 6)
        self.assertTrue(all(str(n) in ['TYDM', 'RTIC', 'TYRT', 'VRYV', 'HIVR', 'YISY'] for n in rf))

        print('>> Passed test_reading_frames()')

    def test_all_orfs(self):
        dna = Dna("ATGAAATTATGAATGAGCCTCAGCTGAAGCATCGCGCATCAGACTACGCTCAGACTCAGACTCAGCATTATAGTGAATGTTAATAAATAAAATAA")
        dna.read_genetic_code('src/tests/files/genetic_code.txt')

        orfs = dna.all_orfs()
        self.assertEqual(orfs[0].get_seq(), 'Protein: MKL')
        self.assertEqual(orfs[1].get_seq(), 'Protein: MSLS')
        self.assertEqual(orfs[2].get_seq(), 'Protein: MLSLSLSVV')
        self.assertEqual(orfs[3].get_seq(), 'Protein: MRDASAEAHS')
        self.assertEqual(orfs[4].get_seq(), 'Protein: MNEPQLKHRASDYAQTQTQHYSEC')
        
        for i in range(0, len(orfs) - 1):
            self.assertTrue(len(orfs[i]) <= len(orfs[i+1]))

        orfs_ten = dna.all_orfs(10)
        self.assertEqual(orfs_ten[0].get_seq(), 'Protein: MRDASAEAHS')
        self.assertEqual(orfs_ten[1].get_seq(), 'Protein: MNEPQLKHRASDYAQTQTQHYSEC')

        print('>> Passed test_all_orfs()')


if __name__ == '__main__':
    unittest.main()
