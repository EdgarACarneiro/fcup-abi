import unittest

from seq_align import MyBlast


class test_MyBlast(unittest.TestCase):

    blast = MyBlast("../testfiles/seqBlast.txt", 11)

    query1 = MyBlast.read_query("../testfiles/query1.fasta")
    query2 = MyBlast.read_query("../testfiles/query2.fasta")

    def test_read_database(self):
        db = MyBlast.read_database("../testfiles/seqBlast.txt")

        self.assertEqual(db[3], 'atagctcgatgcttagatctcgcgtatgctgctagataagagctgctgagctgatcggatgcctcgcgctcgcgcgctgaggctcggatagctagctgagcgctcgatagcgcgttcgctggatcgcgtatagcgctgaagctcccggctagctgtctgtaaatcggatctgatctcgctctatact')

    def test_read_query(self):
        query = MyBlast.read_query("../testfiles/query2.fasta")

        self.assertEqual(query, 'cgacgacgacgacgaatgatg')

    def test_add_seq_to_db(self):
        blast = MyBlast(None, 11)

        self.assertEqual(0, len(blast.database))
        blast.add_sequence_database('test sequence')
        self.assertEqual(1, len(blast.database))

    def test_build_map(self):
        self.blast.build_map(self.query2)

        self.assertEqual(self.blast.mapping['cgacgacgacg'], [0, 3])
        self.assertEqual(self.blast.mapping['gacgacgaatg'], [7])

    def test_get_hits(self):
        self.blast.build_map(self.query2)
        hits = self.blast.get_hits(self.blast.database[4])

        self.assertEqual(hits[0], (0, 0))
        self.assertEqual(hits[3], (4, 1))
        self.assertEqual(hits[4], (2, 2))
        self.assertEqual(hits[7], (1, 4))
        self.assertEqual(hits[10], (6, 6))

    def test_extends_hit(self):
        self.blast.build_map(self.query2)
        hits = self.blast.get_hits(self.blast.database[4])

        self.assertEqual(
            self.blast.extends_hit(
                self.blast.database[4], hits[10], self.query2), # hit = (6,6)
            (0, 0, 21, 21))
        self.assertEqual(
            self.blast.extends_hit(
                self.blast.database[4], hits[0], self.query2), # hit = (0,0)
            (0, 0, 21, 21))

    def test_hit_best_score(self):
        self.blast.build_map(self.query2)
        bs = self.blast.hit_best_score(self.blast.database[4], self.query2)

        self.assertEqual(bs, (0, 0, 21, 21))

    def test_best_alignment(self):
        self.assertEqual(self.blast.best_alignment(self.query2), (0, 0, 21, 21, 4))
        self.assertEqual(self.blast.best_alignment(self.query1), (1, 38, 149, 108, 3))


if __name__ == '__main__':
    unittest.main()
