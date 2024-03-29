import unittest

from ...bioseq import Protein
from ...seq_align import MyAlign


class test_MyAlign(unittest.TestCase):

    align_1 = ["ATGA-A", "AA-AT-"]
    align_protein = [Protein("VJKK"), Protein("JRSK"), Protein("VRSK")]

    def test_create_my_align(self):
        self.assertRaises(Exception, MyAlign, ["AA", Protein("AA")])

        try:
            MyAlign(self.align_1)
            MyAlign(self.align_protein)
        except Exception:
            self.fail("MyAlign() raised Exception unexpectedly!")

        print('>> Passed test_create_my_align()')

    def test_basics(self):
        align = MyAlign(self.align_1)
        align_p = MyAlign(self.align_protein)

        self.assertEqual(str(align), "\nATGA-A\nAA-AT-")
        self.assertEqual(len(align), 6)
        self.assertEqual(align[0], "ATGA-A")
        self.assertEqual(align[1, 3], "A")
        self.assertEqual(align.num_seqs(), 2)
        self.assertEqual(align.get_seqs(), self.align_1)
        self.assertEqual(align.get_align_type(), 'str')

        self.assertEqual(str(align_p), "\nVJKK\nJRSK\nVRSK")
        self.assertEqual(len(align_p), 4)
        self.assertEqual(align_p[2], "VRSK")
        self.assertEqual(align_p[1, 3], "K")
        self.assertEqual(align_p.num_seqs(), 3)
        self.assertEqual(align_p.get_seqs(),
                         [str(prot) for prot in self.align_protein])
        self.assertEqual(align_p.get_align_type(), 'Protein')

        print('>> Passed MyAlign::test_basics()')

    def test_num_column(self):
        align = MyAlign(self.align_1)
        align_p = MyAlign(self.align_protein)

        self.assertEqual("".join(align.column(1)), "TA")
        self.assertEqual("".join(align_p.column(3)), "KKK")

        print('>> Passed test_num_column()')

    def test_consensus(self):
        align = MyAlign(self.align_1)
        align_p = MyAlign(self.align_protein)

        self.assertEqual(align.consensus(), "ATGATA")
        self.assertEqual(align_p.consensus(), "VRSK")

        print('>> Passed test_consensus()')


if __name__ == '__main__':
    unittest.main()
