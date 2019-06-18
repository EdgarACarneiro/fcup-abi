from bioseq import BioSeq, Seq, Dna, Rna, Protein
from seq_align import MyBlast, MultipleAlignment, SubstMatrix
from phylogenetics import UPGMA
import re


class Pipeline:

    query_seq: Seq
    query_id: str
    database: [(str, Seq)] # Used list instead of dict since order is important

    def __init__(self, input_fasta, database_fasta):
        fasta_dic = BioSeq.read_fasta_file(input_fasta)

        if len(fasta_dic) != 1:
            raise Exception(
                "Input dictionary expected to have only one sequence")

        self.query_seq = Pipeline.infer_type(
            list(fasta_dic.values())[0])
        self.query_id = list(fasta_dic.keys())[0]

        db_fasta = BioSeq.read_fasta_file(database_fasta)
        for seq in db_fasta.keys():
            db_fasta[seq] = Pipeline.infer_type(db_fasta[seq])

        self.database = list(map(
            lambda k: (k, db_fasta[k]), db_fasta.keys()))

    @staticmethod
    def infer_type(seq):
        """Given a seq infer its type from its contents"""
        if all([el in Dna.switcher for el in seq]):
            return Dna(seq)

        elif all([el in Rna.switcher for el in seq]):
            return Rna(seq)

        else:
            return Protein(seq)

    @staticmethod
    def get_specie_from_id(seq_id):
        """Get the specie from the given seq_id"""
        return (re.search(r'.*?\[(.*?)\].*', seq_id)).group(1).replace(' ', '_')

    def get_specie(self):
        """Get the query sequence's specie"""
        return Pipeline.get_specie_from_id(self.query_id)

    def get_specie_from_seq(self, seq):
        """Get the sequence's specie. The sequence must be stored in either the
        database or the query sequence."""
        if seq == self.query_seq:
            return self.get_specie()

        for k, v in self.database:
            if v == seq:
                return Pipeline.get_specie_from_id(k)

    def execute(self):
        """Execute the Pipeline"""
        # Create copy database without similar species
        db_copy = [(_id, seq)
                   for _id, seq in self.database
                   if Pipeline.get_specie_from_id(_id) !=
                   self.get_specie()]

        # Creating Blast and populating its database
        blast = MyBlast()
        for _, seq in db_copy:
            blast.add_sequence_database(str(seq))

        # Top 10 Alignments
        top_alignments = blast.best_alignments(str(self.query_seq), 50)
        print([(self.get_specie_from_seq(db_copy[a[4]][1]), a[3]) for a in top_alignments])

        # Mapping the top alignments into the correspondent sequences
        db_seqs = [val for _, val in db_copy]
        best_seqs = [self.query_seq] + list(map(
            lambda align: db_seqs[align[4]], top_alignments))

        # Substitution Matrix and Gap penalty
        align_data = (SubstMatrix.read_submat_file(
            "tests/files/blosum62.mat"), -1)

        # Multiple Sequence Alignment
        msa = MultipleAlignment(best_seqs, align_data).align_consensus()

        # Printing the Multiple Sequence Alignemnt
        print("Multiple Sequence Alignment Result:\n")
        msa.pretty_print()
        print()

        # Producing the Ultrametric Tree
        upgma = UPGMA(best_seqs, align_data)
        tree = upgma.run()

        # Printing the tree with a mapping for the species
        print("Phylogenetics Tree created:\n")
        tree.print_tree({
            idx: self.get_specie_from_seq(seq)
            for idx, seq in enumerate(best_seqs)
        })
        print()

        # 


if __name__ == '__main__':
    p = Pipeline('tests/files/source.fasta', 'tests/files/seqdump.txt')
    p.execute()
    # print(p.database)
