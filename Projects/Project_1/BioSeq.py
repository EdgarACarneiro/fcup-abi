from enum import Enum

class BioSeq:

    class _SeqType(Enum):
        DNA = "DNA"
        RNA = "RNA"
        PROTEIN = "PROTEIN"

    def __init__(self, seq, seq_type="DNA"):
        self.__seq__ = seq
        self.__seq__type__ = self.SeqType(seq_type)
