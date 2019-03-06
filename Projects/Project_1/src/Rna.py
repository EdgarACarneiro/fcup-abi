from .NucleotideChain import NucleotideChain


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
        """Pretty printing of the class information:
        type and sequence"""
        print("* Bio type:\nRNA")
        super().pretty_print()

    def rev_transcription(self):
        """computes the DNA corresponding to the reverse transcription of the RNA sequence"""
        return Rna(self._seq.replace("U", "T"))
