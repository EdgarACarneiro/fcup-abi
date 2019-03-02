from NucleotideChain import NucleotideChain


class Rna(NucleotideChain):
    def __init__(self, seq):
        super().__init__(seq)

    def validate(self):
        return all(n in ["A", "U", "C", "G"] for n in self._seq)