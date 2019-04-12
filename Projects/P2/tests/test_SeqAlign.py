import unittest

from seq_align import read_submat_file,\
                      global_align_multiple_solutions

class test_SeqAlign(unittest.TestCase):

    # def test_init(self):
    #     # self.assertEqual(Rna("UCGA").get_seq(), "UCGA")
    #     # # Launches Exception since sequence contains invalid chars ('T'/ 'S')
    #     # self.assertRaises(Seq.InvalidSequenceException, Rna, "ACTG")
    #     # self.assertRaises(Seq.InvalidSequenceException, Rna, "ACSG")

    def test_read_submat(self):
        sm = read_submat_file('files/blosum62.mat')

        # Choosing random elements to check if sm is loaded correctly
        self.assertEqual(sm['AA'], 4)
        self.assertEqual(sm['AD'], -2)
        self.assertEqual(sm['AC'], 0)
        self.assertEqual(sm['RQ'], 1)
        self.assertEqual(sm['RM'], -1)
        self.assertEqual(sm['NC'], 3)
        self.assertEqual(sm['NP'], -2)
        self.assertEqual(sm['NW'], -4)
        self.assertEqual(sm['DD'], 6)
        self.assertEqual(sm['HP'], -2)
        self.assertEqual(sm['FF'], 6)
        self.assertEqual(sm['VW'], -3)

if __name__ == '__main__':
    unittest.main()
