import unittest

from ..pipeline import Pipeline


class test_pipeline(unittest.TestCase):

    p = Pipeline('src/tests/files/source.fasta', 'src/tests/files/seqdump.txt', 10)

    def test__init__(self):

        self.assertEqual(self.p.query_seq.__class__.__name__, 'Protein')
        self.assertEqual(len(self.p.database), 20)
        for seq in self.p.database:
            self.assertEqual(seq[1].__class__.__name__, 'Protein')

        print('>> Passed Pipeline::test__init__()')

    def test_infer_type(self):
        self.assertEqual(
            self.p.infer_type('ATCG').__class__.__name__, 'Dna')
        self.assertEqual(
            self.p.infer_type('AUAU').__class__.__name__, 'Rna')
        self.assertEqual(
            self.p.infer_type('KLED').__class__.__name__, 'Protein')

        print('>> Passed test_infer_type()')

    def test_get_specie(self):
        self.assertEqual(self.p.get_specie(), 'Homo_sapiens')
        self.assertEqual(
            self.p.get_specie_from_id(self.p.database[3][0]), 'Gorilla_gorilla_gorilla')
        self.assertEqual(
            self.p.get_specie_from_seq(self.p.query_seq), 'Homo_sapiens')
        self.assertEqual(
            self.p.get_specie_from_seq(self.p.database[3][1]), 'Gorilla_gorilla_gorilla')

        print('>> Passed test_get_specie()')

    def test_change_alignment_settings(self):
        self.assertEqual(self.p.align_config[1], -8)

        self.p.change_alignment_settings(None, -1)

        self.assertEqual(self.p.align_config[1], -1)
        self.assertEqual(self.p.align_config[0], None)

        print('>> Passed test_change_alignment_settings()')


if __name__ == '__main__':
    unittest.main()
