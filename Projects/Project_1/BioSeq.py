from enum import Enum
from abc import ABC, abstractmethod


class BioSeq:
    """Class that represents a Bio Sequence: can either be a Dna, Rna or Protein"""

    class _SeqType(Enum):
        """Possible Bio Sequence types"""
        DNA = "DNA"
        RNA = "RNA"
        PROTEIN = "PROTEIN"

    class BioTypeException(Exception):
        """Error raised when bio sequence type is not recognized"""
        pass

    def __init__(self, seq, seq_type="DNA"):
        seq_type = self._SeqType(seq_type)

        if seq_type is self._SeqType.DNA:
            __bio_seq = _Dna(seq)

        elif seq_type is self._SeqType.RNA:
            __bio_seq = _Rna(seq)

        elif seq_type is self._SeqType.PROTEIN:
            __bio_seq = _Protein(seq)

        else:
            raise self.BioTypeException


class _Seq(ABC):
    """Abstract class that represents a sequence"""

    def __init__(self, seq):
        self.__seq = seq


class _Dna(_Seq):
    switcher = {"A": "T",
                "T": "A",
                "G": "C",
                "C": "G"}

    def __init__(self, seq):
        self.Seq.__init__(self, seq)


class _Rna(_Seq):
    def __init__(self, seq):
        self.Seq.__init__(self, seq)


class _Protein(_Seq):
    def __init__(self, seq):
        self.Seq.__init__(self, seq)
