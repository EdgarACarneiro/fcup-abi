import unittest

from seq_align import read_submat_file,\
                        subst_matrix,\
                        pretty_matrix,\
                        global_align_multiple_solutions

class test_SeqAlign(unittest.TestCase):

    seq1 = 'GATTACA'
    seq2 = 'GCATGCT'

    def test_read_submat(self):
        sm = read_submat_file('tests/files/blosum62.mat')

        # Choosing random elements to check if sm is loaded correctly
        self.assertEqual(sm['AA'], 4)
        self.assertEqual(sm['AD'], -2)
        self.assertEqual(sm['AC'], 0)
        self.assertEqual(sm['RQ'], 1)
        self.assertEqual(sm['RM'], -1)
        self.assertEqual(sm['NN'], 6)
        self.assertEqual(sm['NP'], -2)
        self.assertEqual(sm['NW'], -4)
        self.assertEqual(sm['DD'], 6)
        self.assertEqual(sm['HP'], -2)
        self.assertEqual(sm['YY'], 7)
        self.assertEqual(sm['VW'], -3)

    def test_subst_matrix(self):
        sm = subst_matrix("ACGT", 1, -1)
        self.assertEqual(sm,
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

    def test_global_align_multiple_solutions(self):
        sm = subst_matrix("ACGT", 1, -1)
        ga_score, ga_trace = global_align_multiple_solutions(self.seq1, self.seq2, sm, -3)

        self.assertEqual(ga_score,
            [
                [  0,  -3,  -6, -9, -12, -15, -18, -21],
                [ -3,   1,  -2, -5,  -8, -11, -14, -17],
                [ -6,  -2,   0, -1,  -4,  -7, -10, -13],
                [ -9,  -5,  -3, -1,   0,  -3,  -6,  -9],
                [-12,  -8,  -6, -4,   0,  -1,  -4,  -5],
                [-15, -11,  -9, -5,  -3,  -1,  -2,  -5],
                [-18, -14, -10, -8,  -6,  -4,   0,  -3],
                [-21, -17, -13, -9,  -9,  -7,  -3,  -1]
            ])
        self.assertEqual(ga_trace,
            [
                [[0], [3], [3],    [3],    [3],    [3],    [3],    [3]   ],
                [[2], [0], [2],    [2],    [2],    [0, 2], [2],    [2]   ],
                [[2], [1], [0],    [0],    [2],    [2],    [2],    [2]   ],
                [[2], [1], [0, 1], [0],    [0],    [2],    [2],    [0, 2]],
                [[2], [1], [0, 1], [0, 1], [0],    [0],    [0, 2], [0]   ],
                [[2], [1], [0, 1], [0],    [1],    [0],    [0],    [0, 2]],
                [[2], [1], [0],    [1],    [0, 1], [0, 1], [0],    [0, 2]],
                [[2], [1], [1],    [0],    [0, 1], [0, 1], [1],    [0]   ]
            ])


if __name__ == '__main__':
    unittest.main()
