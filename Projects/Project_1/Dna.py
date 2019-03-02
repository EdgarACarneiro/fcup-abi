from NucleotideChain import NucleotideChain
from Rna import Rna


class Dna(NucleotideChain):

    switcher = {"A": "T",
                "T": "A",
                "G": "C",
                "C": "G"}

    # TODO: Validate in sequences on creation?
    def __init__(self, seq):
        super().__init__(seq)
        self._genetic_code = {}

    def transcription(self):
        return Rna(self._seq.replace("T", "U"))

    def read_genetic_code(self, file_name):
        self._genetic_code = 
            {line[1:4]: line[7] for line in self.readFile(file_name)}
