from .NucleotideChain import NucleotideChain
from .Rna import Rna
from .Protein import Protein


class Dna(NucleotideChain):

    switcher = {"A": "T",
                "T": "A",
                "G": "C",
                "C": "G"}

    def __init__(self, seq):
        super().__init__(seq)
        self._genetic_code = {}

    def __str__(self):
        return 'DNA: ' + self._seq

    def pretty_print(self):
        """Pretty printing of the class information:
        type, sequence and genetic code dictionary"""
        print("* Bio type:\nDNA")
        super().pretty_print()

        print("* Genetic Code Dictionary:")
        for k, v in self._genetic_code.items():
            print(k + ' - ' + v)

    def transcription(self):
        """computes the RNA corresponding to the transcription of the DNA sequence"""
        return Rna(self._seq.replace("T", "U"))

    def read_genetic_code(self, file_name):
        """Reads and stores the genetic code dicitonary from the given file"""
        fd = self.read_file(file_name)
        self._genetic_code = {line[1:4]: line[7] for line in fd}
        fd.close()

    def get_genetic_code(self):
        return self._genetic_code

    def set_genetic_code(self, gc):
        self._genetic_code = gc

    def translate(self, iniPos=0):
        """Translate the stored dna sequence using the stored genetic code
        and returns the resultant Protein."""
        trans = ""

        for i in range(iniPos, len(self._seq) - 2, 3):
            trans += self._genetic_code[self._seq[i: i + 3]]
        return Protein(trans)

    def codon_usage(self, aa, iniPos=0):
        """Provides the frequency of each codon encoding a given aminoacid,
        in the stored dna sequence with the stored dictionary"""
        freq = {}
        total = 0

        for k, v in self._genetic_code.items():
            if v is aa:
                freq[k] = 0

        for i in range(iniPos, len(self._seq) - 2, 3):
            if self._seq[i: i + 3] in freq:
                freq[self._seq[i: i + 3]] += 1
                total += 1

        if total > 0:
            # Filter the dictionary
            return {k: v/total for (k, v) in freq.items() if v > 0}
        else:
            return freq

    def reading_frames(self):
        """Compute all possible reading frames of the stored dna sequence,
        using the stored genetic code (includes the reverse complement)"""
        rc = Dna(self.reverse_complement())
        rc.set_genetic_code(self.get_genetic_code())
        return [self.translate(i) for i in range(0, 3)] + [rc.translate(i) for i in range(0, 3)]

    def __all_orfs_unordered(self):
        """Computes all possible proteins for all open reading frames and
        returnd and unordered list"""
        return [p for rf in self.reading_frames() for p in rf.all_proteins_rf()]

    def all_orfs(self, minsize=0):
        """Computes all possible proteins for all open reading frames.
        Returns ordered list of proteins with minimum size"""
        return sorted([Protein(el) for el in self.__all_orfs_unordered() if len(el) >= minsize], key=lambda prot: len(prot))
