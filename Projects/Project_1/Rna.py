from Seq import Seq


class Rna(Seq):
    def __init__(self, seq):
        super().__init__(seq)

    def validate(self):
        return self._seq in ["A", "U", "C", "G"]