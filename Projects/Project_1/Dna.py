from Seq import Seq


class Dna(Seq):

    def __init__(self, seq):
        super().__init__(seq)
        self.switcher = {"A": "T",
                         "T": "A",
                         "G": "C",
                         "C": "G"}

    def validate(self):
        return self._seq in self.switcher