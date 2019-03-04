from NucleotideChain import NucleotideChain


class Rna(NucleotideChain):

    switcher = {"A": "U",
                "U": "A",
                "G": "C",
                "C": "G"}

    def __init__(self, seq):
        super().__init__(seq)

    def __str__(self):
        return 'RNA: ' + self._seq

    def pretty_print(self):
        print("* Bio type:\nRNA")
        super().pretty_print()    
