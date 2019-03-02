from abc import ABC
from Seq import Seq


class NucleotideChain(Seq, ABC):

    def __init__(self, seq):
        super().__init__(seq)

    def __gc_percent_helper(seq):
        """"Return the percentage of G and C nucleotides in the given sequence"""
        return sum(1 for val in dna_seq if val is "G" or val is "C") / len(seq)

    def gc_percent(self):
        """Return the percentage of G and C nucleotides in the stored Bio Sequence.
        Genes are tipically found in GC-rich regions of the genome"""
        return self.__gc_percent_helper(self._seq)

    def gc_percent_sub_seq(k):
        """Return the percentage of G and C nucleotides in subsequences of size k
        belonging to the stored bio sequence"""
        return list(map(lambda ss: self.__gc_percent_helper(ss),
                        [self._seq[i: i+k] for i in range(0, len(self._seq), k)]))
