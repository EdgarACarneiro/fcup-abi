from .Seq import Seq


class Protein(Seq):

    def __init__(self, seq):
        super().__init__(seq)

    def validate(self):
        return all(n.isalpha() or n is '_' for n in self._seq)

    def __str__(self):
        return 'Protein: ' + self._seq

    def pretty_print(self):
        """Pretty printing of the class information:
        type and sequence"""
        print("* Bio type:\nProtein")
        super().pretty_print()

    def all_proteins_rf(self):
        """Computes all possible proteins in the stored aminoacid sequence
        Complexity: O(log(n))"""
        proteins = []
        curr = []
        begin = False

        for aa in self._seq:
            if aa is "M":
                curr.append("M")
                begin = True
                continue

            if aa is "_" and begin:
                for i in range(0, len(curr)):
                    seq = ""
                    for j in range(i, len(curr)):
                        seq += curr[j]
                    proteins.append(seq)
                curr = []
                begin = False

            if begin:
                curr[len(curr) - 1] += aa

        return proteins
