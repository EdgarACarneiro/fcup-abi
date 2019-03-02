from abc import ABC, abstractmethod


class Seq(ABC):
    """Abstract class that represents a sequence"""

    class InvalidSequenceException(Exception):
        """Error raised when sequence contains invalid characters"""
        pass

    def __init__(self, seq):
        self._seq = seq

    @abstractmethod
    def validate(self):
        pass

    def freq(self):
        """Computes the frequency of the symbols in the stored sequence and 
        returns a sorted dictionary"""
        symbols = {}
        for base in self._seq:
            if (base in symbols):
                symbols[base] += 1
            else:
                symbols[base] = 0

        return sorted(symbols.items(),
                      key=lambda symbol: symbol[1],
                      reverse=True)

    def gcPercent(self):
        """Return the percentage og G and C nucleotides in a Bio Sequence.
        Genes are tipically found in GC-rich regions of the genome"""
        return sum(1 for val in dna_seq if val is "G" or val is "C") / len(dna_seq)