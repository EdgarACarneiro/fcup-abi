from enum import Enum

from Dna import Dna
from Rna import Rna
from Protein import Protein


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

        bio_type = self._BioType(seq_type.upper())

        if bio_type is self._BioType.DNA:
            __bio_seq = Dna(seq)
        elif seq_type is self._BioType.RNA:
            __bio_seq = Rna(seq)
        else:
            __bio_seq = Protein(seq)

    def validate_seq(self):
        return __bio_seq.validate()


def main():
    test = BioSeq("ATCACA", "dna") 
    assert test.validate_seq(), "Invalid Sequence"

if __name__ == "__main__":
    main()