from enum import Enum

from Dna import Dna
from Rna import Rna
from Protein import Protein


class BioSeq:
    """Class that represents a Bio Sequence: can either be a Dna,
    Rna or Protein"""

    class BioType(Enum):
        """Possible Bio Sequence types"""
        DNA = "DNA"
        RNA = "RNA"
        PROTEIN = "PROTEIN"

    class BioTypeException(Exception):
        """Error raised when bio sequence type is not recognized"""
        pass

    @staticmethod
    def createBioSeq(seq, seq_type="DNA"):
        """Create a BioSeq instance according to the given type and
        returns it"""

        # Assuring it is either Dna, Rna or a Protein
        if seq_type.upper() not in [el.value for el in list(BioSeq.BioType)]:
            raise self.BioTypeException

        bio_type = BioSeq.BioType(seq_type.upper())

        if bio_type is BioSeq.BioType.DNA:
            return Dna(seq)
        elif seq_type is BioSeq.BioType.RNA:
            return Rna(seq)
        return Protein(seq)


def main():
    dna = BioSeq.createBioSeq("AtcACA", "dna")
    assert dna.validate(), "Invalid Sequence"


if __name__ == "__main__":
    main()