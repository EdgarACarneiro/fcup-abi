from Seq import Seq


class Dna(Seq):
    switcher = {"A": "T",
                "T": "A",
                "G": "C",
                "C": "G"}

    def __init__(self, seq):
        super().__init__(seq)

    def validate(self):
        return self.__seq in switcher