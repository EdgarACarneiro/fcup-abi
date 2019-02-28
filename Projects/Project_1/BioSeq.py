from enum import Enum

import Dna
import Rna
import Protein


class BioSeq:
    """Class that represents a Bio Sequence: can either be a Dna,
    Rna or Protein"""

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