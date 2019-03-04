from NucleotideChain import NucleotideChain
from Rna import Rna
from Protein import Protein


class Dna(NucleotideChain):

    switcher = {"A": "T",
                "T": "A",
                "G": "C",
                "C": "G"}

    # TODO: Validate in sequences on creation?
    def __init__(self, seq):
        super().__init__(seq)
        self._genetic_code = {}

    def transcription(self):
        return Rna(self._seq.replace("T", "U"))

    def read_genetic_code(self, file_name):
        self._genetic_code = {line[1:4]: line[7] for line in self.readFile(file_name)}

    def translate(iniPos = 0):
        """Translate the stored dna sequence using the stored dictionary
        and returns the resultant Protein"""
        trans = ""
        
        for i in range(iniPos, len(seq) - 2, 3):
            trans += dictionary[seq[i: i + 3]]
        return Protein(trans)

    def codonUsage(aa, iniPos = 0):
        """Provides the frequency of each codon encoding a given aminoacid,
        in the stored dna sequence with the stored dictionary"""
        freq = {}
        total = 0
        
        for k,v in self._genetic_code.items():
            if v is aa:
                freq[k] = 0
        
        for i in range(iniPos, len(seq) - 2, 3):
            if self._seq[i: i + 3] in freq:
                freq[self._seq[i: i + 3]] += 1
                total += 1
                
        if total > 0:
            return {k: v/total for (k,v) in freq.items() if v > 0} # Filter the dictionary
        else:
            return freq