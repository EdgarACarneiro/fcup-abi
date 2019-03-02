from NucleotideChain import NucleotideChain


class Rna(NucleotideChain):

    switcher = {"A": "U",
                "U": "A",
                "G": "C",
                "C": "G"}

    def __init__(self, seq):
        super().__init__(seq)
