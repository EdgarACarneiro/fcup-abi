from enum import Enum
from abc import ABC, abstractmethod


class BioSeq:
    """Class that represents a Bio Sequence: can either be a Dna, Rna or Protein"""

    class _BioType(Enum):
        """Possible Bio Sequence types"""
        DNA = "DNA"
        RNA = "RNA"
        PROTEIN = "PROTEIN"

    class BioTypeException(Exception):
        """Error raised when bio sequence type is not recognized"""
        pass

    def __init__(self, seq, seq_type="DNA"):
        # Assuring it is either Dna, Rna or a Protein
        if seq_type.upper() not in [el.value for el in list(self._BioType)]:
            raise self.BioTypeException

        bio_type = self._BioType(seq_type)

        if bio_type is self._BioType.DNA:
            __bio_seq = _Dna(seq)
        elif seq_type is self._BioType.RNA:
            __bio_seq = _Rna(seq)
        else:
            __bio_seq = _Protein(seq)

# Não preciso disto -> podem herdar de bioseq fds
class _Seq(ABC):
    """Abstract class that represents a sequence"""

    class InvalidSequenceException(Exception):
        """Error raised when sequence contains invalid characters"""
        pass

    def __init__(self, seq):
        self.__seq = seq


class _Dna(_Seq):
    switcher = {"A": "T",
                "T": "A",
                "G": "C",
                "C": "G"}

    def __init__(self, seq):
        super().__init__(seq)


class _Rna(_Seq):
    def __init__(self, seq):
        super().__init__(seq)


class _Protein(_Seq):
    def __init__(self, seq):
        super().__init__(seq)
