from .NucleotideChain import NucleotideChain


class Rna(NucleotideChain):

    switcher = {"A": "U",
                "U": "A",
                "G": "C",
                "C": "G"}

    def __init__(self, seq):
        super().__init__(seq)

    def get_seq(self):
        return 'RNA: ' + self._seq

    def pretty_print(self):
        """Pretty printing of the class information:
        type and sequence"""
        print("* Bio type:\nRNA")
        super().pretty_print()

    def rev_transcription(self):
        """computes the DNA corresponding to the reverse transcription of the RNA sequence"""
        from .Dna import Dna  # Import done here because of circular dependencies
        return Dna(self._seq.replace("U", "T"))
