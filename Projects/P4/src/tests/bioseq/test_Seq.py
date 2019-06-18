import unittest

from ...bioseq import Seq
from ...bioseq import BioSeq


class test_Seq(unittest.TestCase):

    def test_init(self):
        # Launches Exception since class is abstract
        self.assertRaises(Exception, Seq, "ACTG")

        print('>> Passed Seq::test_init()')

    def test_basics(self):
        seq = BioSeq.create_bio_seq("ACTGTCATGAT")
        self.assertEqual(str(seq), "ACTGTCATGAT")
        self.assertEqual(len(seq), len("ACTGTCATGAT"))

        print('>> Passed Seq::test_basics()')

    def test_freq(self):
        seq = BioSeq.create_bio_seq("ACTGTCATAT")
        self.assertEqual(seq.freq_symbols(), {"T": 4, "A": 3, "C": 2, "G": 1})

        print('>> Passed test_freq()')

    def test_read_write(self):
        seq = BioSeq.create_bio_seq("ACTGTCATAT")
        len_bwr = len(seq)

        seq.write_sequence('src/tests/files/test_write.csv')
        seq.read_sequence('src/tests/files/seq.csv')
        self.assertGreater(len(seq), len_bwr)
        self.assertEqual(2591, len(seq))

        seq.read_sequence('src/tests/files/test_write.csv')
        self.assertEqual(len_bwr, len(seq))

        print('>> Passed test_read_write()')

    def test_save_load(self):
        seq = BioSeq.create_bio_seq("ACTGTCATAT")
        seq.save('src/tests/files/test_save_load.csv')
        loaded_seq = BioSeq.load('src/tests/files/test_save_load.csv')
        self.assertEqual(seq.get_seq(), loaded_seq.get_seq())
        self.assertEqual(seq.get_genetic_code(), loaded_seq.get_genetic_code())
        self.assertEqual(seq.__class__.__name__, loaded_seq.__class__.__name__)

        print('>> Passed test_save_load()')


if __name__ == '__main__':
    unittest.main()
