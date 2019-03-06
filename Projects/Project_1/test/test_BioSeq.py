import unittest

from src import BioSeq
from src import Dna


class test_BioSeq(unittest.TestCase):

    def test_create_bio_seq(self):
        self.assertIsInstance(BioSeq.create_bio_seq("ATCGTAT"), Dna)

if __name__ == '__main__':
    unittest.main()