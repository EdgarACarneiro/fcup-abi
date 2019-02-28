class _Dna(_Seq):
    switcher = {"A": "T",
                "T": "A",
                "G": "C",
                "C": "G"}

    def __init__(self, seq):
        super().__init__(seq)
