import unittest

from ...seq_align import SubstMatrix

class test_SubstMatrix(unittest.TestCase):

    sm = SubstMatrix.read_submat_file('src/tests/files/blosum62.mat')

    def test_read_submat(self):
        # Choosing random elements to check if sm is loaded correctly
        self.assertEqual(self.sm['AA'], 4)
        self.assertEqual(self.sm['AD'], -2)
        self.assertEqual(self.sm['AC'], 0)
        self.assertEqual(self.sm['RQ'], 1)
        self.assertEqual(self.sm['RM'], -1)
        self.assertEqual(self.sm['NN'], 6)
        self.assertEqual(self.sm['NP'], -2)
        self.assertEqual(self.sm['NW'], -4)
        self.assertEqual(self.sm['DD'], 6)
        self.assertEqual(self.sm['HP'], -2)
        self.assertEqual(self.sm['YY'], 7)
        self.assertEqual(self.sm['VW'], -3)

        print('>> Passed test_read_submat()')

    def test_create_submat(self):
        sm = SubstMatrix.create_submat("ACGT", 1, -1)
        self.assertEqual(sm.sm,
                         {'AA': 1,
                          'AC': -1,
                          'AG': -1,
                          'AT': -1,
                          'CA': -1,
                          'CC': 1,
                          'CG': -1,
                          'CT': -1,
                          'GA': -1,
                          'GC': -1,
                          'GG': 1,
                          'GT': -1,
                          'TA': -1,
                          'TC': -1,
                          'TG': -1,
                          'TT': 1})

        print('>> Passed test_create_submat()')

if __name__ == '__main__':
    unittest.main()
