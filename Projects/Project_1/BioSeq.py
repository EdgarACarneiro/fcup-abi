from enum import Enum
import pickle

from Dna import Dna
from Rna import Rna
from Protein import Protein


class BioSeq:
    """Class that represents a Bio Sequence: can either be a Dna,
    Rna or Protein"""

    class BioType(Enum):
        """Possible Bio Sequence types"""
        DNA = "DNA"
        RNA = "RNA"
        PROTEIN = "PROTEIN"

    class BioTypeException(Exception):
        """Error raised when bio sequence type is not recognized"""
        pass

    @staticmethod
    def createBioSeq(seq, seq_type="DNA"):
        """Create a BioSeq instance according to the given type and
        returns it"""

        # Assuring it is either Dna, Rna or a Protein
        if seq_type.upper() not in [el.value for el in list(BioSeq.BioType)]:
            raise self.BioTypeException

        bio_type = BioSeq.BioType(seq_type.upper())

        if bio_type is BioSeq.BioType.DNA:
            return Dna(seq)
        elif seq_type is BioSeq.BioType.RNA:
            return Rna(seq)
        return Protein(seq)

    @staticmethod
    def readFastaFile(fileName):
        """Reads a fasta file and returns a dicitonary with the correspondent
        identifier and bio sequence."""
        fasta = {}
        helperSeq = ""
        currFasta = ""
        
        for line in readFile(fileName):
            l = line.strip()
            
            if l is "" and currFasta is not "":
                fasta[currFasta] = helperSeq
                helperSeq = ""
                currFasta = ""
            
            if currFasta is not "":
                helperSeq += l
            
            if (l.strip()[0:1] is ">"):
                currFasta = line[1:len(l)]
                
        return fasta

    @staticmethod
    def load(fileName):
        return pickle.load(open(fileName, mode='rb'))

def main():
    dna = BioSeq.createBioSeq("ATGAAATTATGAATGAGCCTCAGCTGAAGCATCGCGCATCAGACTACGCTCAGACTCAGACTCAGCATTATAGTGAATGTTAATAAATAAAATAA", "dna")
    assert dna.validate(), "Invalid Sequence"
    print(dna.reverse_complement())
    print(dna.transcription().reverse_complement())
    dna.read_genetic_code('../../Classes/files/genetic_code.txt')
    for p in dna.reading_frames():
        print(p.get_seq())
    print(dna.all_orfs(10))
    dna.pretty_print()
    # print(dna.load('test.csv'))


if __name__ == "__main__":
    main()