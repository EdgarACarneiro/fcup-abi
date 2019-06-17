from enum import Enum
import pickle

from .Dna import Dna
from .Rna import Rna
from .Protein import Protein


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
    def create_bio_seq(seq, seq_type="DNA"):
        """Create a BioSeq instance according to the given type and
        returns it"""

        # Assuring it is either Dna, Rna or a Protein
        if seq_type.upper() not in [el.value for el in list(BioSeq.BioType)]:
            raise BioSeq.BioTypeException

        bio_type = BioSeq.BioType(seq_type.upper())

        if bio_type is BioSeq.BioType.DNA:
            return Dna(seq)
        elif bio_type is BioSeq.BioType.RNA:
            return Rna(seq)
        return Protein(seq)

    @staticmethod
    def read_fasta_file(file_name):
        """Reads a fasta file and returns a dicitonary with the correspondent
        identifier and bio sequence."""
        fasta = {}
        helper_seq = ""
        curr_fasta = ""
        fd = open(file_name, "r")

        for line in fd:
            l = line.strip()

            if (l is "" or l.strip()[0:1] is ">") and curr_fasta is not "":
                fasta[curr_fasta] = helper_seq
                helper_seq = ""
                curr_fasta = ""

            if curr_fasta is not "":
                helper_seq += l

            if (l.strip()[0:1] is ">"):
                curr_fasta = line[1:len(l)]

        fd.close()
        return fasta

    @staticmethod
    def load(file_name):
        """Loads a bio sequence from a given file and returns it"""
        fd = open(file_name, mode='rb')
        load = pickle.load(fd)
        fd.close()
        return load
