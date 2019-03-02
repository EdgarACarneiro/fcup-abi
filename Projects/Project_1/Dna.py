from NucleotideChain import NucleotideChain
from Rna import Rna


class Dna(NucleotideChain):

    # TODO: Validate in sequences on creation?
    def __init__(self, seq):
        super().__init__(seq)
        self.switcher = {"A": "T",
                         "T": "A",
                         "G": "C",
                         "C": "G"}

    def validate(self):
        return all(n in self.switcher for n in self._seq)

    def transcription(self):
        return Rna(self._seq.replace("T", "U"))