import unittest

from bioseq import BioSeq
from seq_align import read_submat_file,\
                        subst_matrix,\
                        pretty_matrix,\
                        global_align_multiple_solutions,\
                        recover_global_align_multiple_solutions,\
                        local_align_multiple_solutions,\
                        recover_local_align_multiple_solutions,\
                        compare_pairwise_global_align,\
                        compare_pairwise_local_align,\
                        compare_pairwise_num_global_align,\
                        compare_pairwise_num_local_align

class test_SeqAlign(unittest.TestCase):

    sm = read_submat_file('tests/files/blosum62.mat')
    sm_dna = subst_matrix("ACGT", 1, -1)

    seq1 = 'GATTACA'
    seq2 = 'GCATGCT'

    slides_seq1 = 'HGWAG'
    slides_seq2 = 'PHSWG'

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

        print('>> Passed test_subst_matrix()')

    def test_global_align_multiple_solutions(self):
        ga_score, ga_trace = global_align_multiple_solutions(self.seq1, self.seq2, self.sm_dna, -3)

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
                [[2], [1], [3],    [3],    [3],    [1, 3], [3],    [3]   ],
                [[2], [2], [1],    [1],    [3],    [3],    [3],    [3]   ],
                [[2], [2], [1, 2], [1],    [1],    [3],    [3],    [1, 3]],
                [[2], [2], [1, 2], [1, 2], [1],    [1],    [1, 3], [1]   ],
                [[2], [2], [1, 2], [1],    [2],    [1],    [1],    [1, 3]],
                [[2], [2], [1],    [2],    [1, 2], [1, 2], [1],    [1, 3]],
                [[2], [2], [2],    [1],    [1, 2], [1, 2], [2],    [1]   ]
            ])

        print('>> Passed test_global_align_multiple_solutions()')

    def test_recover_global_align_multiple_solutions(self):
        _, ga_trace = global_align_multiple_solutions(self.slides_seq1, self.slides_seq2, self.sm, -3)
        rga = recover_global_align_multiple_solutions(ga_trace, self.slides_seq1, self.slides_seq2)
        seq1_alignments = [align[0] for align in rga]
        seq2_alignments = [align[1] for align in rga]

        # Classes Example
        self.assertTrue('-HGWAG' in seq1_alignments)
        self.assertTrue('PHSW-G' in seq2_alignments)

        # C2 example
        _, ga_trace = global_align_multiple_solutions(self.seq1, self.seq2, self.sm_dna, -1)
        rga = recover_global_align_multiple_solutions(ga_trace, self.seq1, self.seq2)
        seq1_alignments = [align[0] for align in rga]
        seq2_alignments = [align[1] for align in rga]

        self.assertEqual(len(seq1_alignments), 3)
        self.assertTrue('G-ATTACA' in seq1_alignments)
        self.assertTrue('GCA-TGCT' in seq2_alignments)
        self.assertTrue('GCAT-GCT' in seq2_alignments)
        self.assertTrue('GCATG-CT' in seq2_alignments)

        seqs = BioSeq.read_fasta_file('tests/files/protein_sequences.fas')
        _, ga_trace = global_align_multiple_solutions(seqs['sp|C1F111'], seqs['sp|B7JC18'], self.sm, -3)
        rga = recover_global_align_multiple_solutions(ga_trace, seqs['sp|C1F111'], seqs['sp|B7JC18'])

        # 5760 global optimal alignments between sp|C1F111 & sp|B7JC18
        self.assertEqual(len(rga), 5760)

        print('>> Passed test_recover_global_align_multiple_solutions()')

    def test_local_align_multiple_solutions(self):
        # Example presented in Class slides
        ga_score, ga_trace, max_score = local_align_multiple_solutions(self.slides_seq2, self.slides_seq1, self.sm, -8)
        self.assertEqual(ga_score,
            [
                [0, 0, 0,  0,  0,  0],
                [0, 0, 0,  0,  0,  0],
                [0, 8, 0,  0,  0,  0],
                [0, 0, 8,  0,  1,  0],
                [0, 0, 0, 19, 11,  3],
                [0, 0, 6, 11, 19, 17]
            ])
        self.assertEqual(ga_trace,
            [
                [[0], [0], [0], [0], [0], [0]],
                [[0], [0], [0], [0], [0], [0]],
                [[0], [1], [0], [0], [0], [0]],
                [[0], [0], [1], [0], [1], [0]],
                [[0], [0], [0], [1], [3], [3]],
                [[0], [0], [1], [2], [1], [1]]
            ])
        self.assertEqual(max_score, 19)

        # Example presented in Class slides
        _, ga_trace, max_score = local_align_multiple_solutions(self.seq1, self.seq2, self.sm_dna, -1)
        self.assertEqual(ga_trace,
            [
                [[0], [0], [0], [0], [0],    [0], [0], [0]],
                [[0], [1], [0], [0], [0],    [1], [0], [0]],
                [[0], [0], [0], [1], [0],    [0], [0], [0]], 
                [[0], [0], [0], [0], [1],    [3], [0], [1]],
                [[0], [0], [0], [0], [1, 2], [1], [0], [1]],
                [[0], [0], [0], [1], [0],    [0], [0], [0]],
                [[0], [0], [1], [0], [0],    [0], [1], [0]],
                [[0], [0], [0], [1], [3],    [0], [0], [0]]
            ])
        self.assertEqual(max_score, 2)

        print('>> Passed test_local_align_multiple_solutions()')

    def test_recover_local_align_multiple_solutions(self):
        # Classes example
        ga_score, ga_trace, _ = local_align_multiple_solutions(self.slides_seq2, self.slides_seq1, self.sm, -8)
        rga = recover_local_align_multiple_solutions(ga_score, ga_trace, self.slides_seq2, self.slides_seq1)

        self.assertEqual(rga, [['HSW', 'HGW'], ['HSWG', 'HGWA']])

        ga_score, ga_trace, _ = local_align_multiple_solutions(self.seq1, self.seq2, self.sm_dna, -1)
        rga = recover_local_align_multiple_solutions(ga_score, ga_trace, self.seq1, self.seq2)

        self.assertEqual(rga, [['AT', 'AT'], ['CA', 'CA']])

        seqs = BioSeq.read_fasta_file('tests/files/protein_sequences.fas')
        ga_score, ga_trace, _ = local_align_multiple_solutions(seqs['sp|C1F111'], seqs['sp|B7JC18'], self.sm, -3)
        rga = recover_local_align_multiple_solutions(ga_score, ga_trace, seqs['sp|C1F111'], seqs['sp|B7JC18'])

        # 4 local optimal alignments between sp|C1F111 & sp|B7JC18
        self.assertEqual(len(rga), 4)

        print('>> Passed test_recover_local_align_multiple_solutions()')

    def test_compare_pairwise_global_align(self):
        seqs = list(BioSeq.read_fasta_file('tests/files/protein_sequences.fas').values())
        cga = compare_pairwise_global_align(seqs, self.sm, -3)

        # Some random values
        self.assertEqual(cga[4][3], 196) # Between sp|C1F111 & sp|B7JC18 - matches previous test
        self.assertEqual(cga[0][0], 1018)
        self.assertEqual(cga[9][6], 514)
        self.assertEqual(cga[8][10], 524)

        print('>> Passed test_compare_pairwise_global_align()')

    def test_compare_pairwise_local_align(self):
        seqs = list(BioSeq.read_fasta_file('tests/files/protein_sequences.fas').values())
        cla = compare_pairwise_local_align(seqs, self.sm, -3)

        # Some random values
        self.assertEqual(cla[4][3], 353) # Between sp|C1F111 & sp|B7JC18 - matches previous test
        self.assertEqual(cla[0][0], 1018)
        self.assertEqual(cla[9][6], 517)
        self.assertEqual(cla[8][10], 524)

        print('>> Passed test_compare_pairwise_local_align()')

    def test_compare_pairwise_num_global_align(self):
        seqs = list(BioSeq.read_fasta_file('tests/files/protein_sequences.fas').values())
        cga = compare_pairwise_num_global_align(seqs, self.sm, -3)

        # Some random values
        self.assertEqual(cga[4][3], 5760) # Between sp|C1F111 & sp|B7JC18 - matches previous test
        self.assertEqual(cga[0][0], 1)
        self.assertEqual(cga[9][6], 288)
        self.assertEqual(cga[8][10], 1728)

        print('>> Passed test_compare_num_pairwise_global_align()')

    def test_compare_pairwise_num_local_align(self):
        seqs = list(BioSeq.read_fasta_file('tests/files/protein_sequences.fas').values())
        cla = compare_pairwise_num_local_align(seqs, self.sm, -3)

        # Some random values
        self.assertEqual(cla[4][3], 4) # Between sp|C1F111 & sp|B7JC18 - matches previous test
        self.assertEqual(cla[0][0], 1)
        self.assertEqual(cla[9][6], 144)
        self.assertEqual(cla[8][10], 1152)

        print('>> Passed test_compare_num_pairwise_local_align()')

if __name__ == '__main__':
    unittest.main()
