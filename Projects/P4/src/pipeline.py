from .bioseq import BioSeq, Seq, Dna, Rna, Protein
from .seq_align import MyBlast, MultipleAlignment, SubstMatrix
from .phylogenetics import UPGMA, MyGraph
import re


class Pipeline:
    """Class responsible for implementing the Pipeline of automatic sequence analysis.´
    First setup the Class using the constructor its methods or by changing the attributes.
    Then, run the pipeline with the execute method."""

    query_seq: Seq
    query_id: str
    # Used list instead of dict since order is important
    database: [(str, Seq)]
    cut: [float]
    align_config: (SubstMatrix, float)

    # Pipeline changeable configs
    TOP = 10
    WORD_SIZE = 3

    def __init__(self, input_fasta, database_fasta, cut, sub_mat_file=None, gap=-8):
        """input_fasta: file with the query fasta.
        database_fasta: file with the fastas that will be the database.
        cut: maximum distance accepted in graph creation"""
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

        if sub_mat_file != None:
            self.align_config =\
                (SubstMatrix.read_submat_file(sub_mat_file), gap)
        else:
            self.align_config = (None, gap)

        self.cut = cut

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

    def change_alignment_settings(self, subst_mat, gap):
        """Change the settings used in the alignment:
        (SubsMatrix and gap penalty)"""
        self.align_config = (subst_mat, gap)

    def execute(self):
        """Execute the Pipeline"""
        print("\n\t:::Step 1 - Create copy database without similar specie sequences:::\n")
        # Create copy database without similar species
        db_copy = [(_id, seq)
                   for _id, seq in self.database
                   if Pipeline.get_specie_from_id(_id) !=
                   self.get_specie()]

        # Print Copy Database Info
        print("Correspondent species of sequences in copy Database:")
        for _id, _ in db_copy:
            print("\t%s" % Pipeline.get_specie_from_id(_id))

        print("\n\t:::Step 2 - Running BLAST and getting top alignments:::\n")
        # Creating Blast and populating its database
        blast = MyBlast(word_size=self.WORD_SIZE)
        for _, seq in db_copy:
            blast.add_sequence_database(str(seq))

        # Top 10 Alignments
        top_alignments = blast.best_alignments(str(self.query_seq), self.TOP)

        # Printing Top Alignments
        print("Top %d Alignments obtained from Blast" % self.TOP)
        for al in top_alignments:
            print(al)

        # Mapping the top alignments into the correspondent sequences
        db_seqs = [val for _, val in db_copy]
        best_seqs = [self.query_seq] + list(map(
            lambda align: db_seqs[align[4]], top_alignments))

        print("\n\t:::Step 3 - Running MSA with the respective top alignments:::\n")
        # Multiple Sequence Alignment
        msa = MultipleAlignment(best_seqs, self.align_config).align_consensus()

        # Printing the Multiple Sequence Alignemnt
        print("Multiple Sequence Alignment Result:\n")
        msa.pretty_print()
        print()

        print("\n\t:::Step 4 - Obtaining the Ultrametric Tree from the top alignments, using UPGMA:::\n")
        # Setting up UPGMA
        upgma = UPGMA(best_seqs, self.align_config)

        # Printing UPGMA distance matrix
        print("Distances Matrix obtained by the UPGMA method:")
        upgma.dists_mat.print_mat()

        # Producing the Ultrametric Tree
        tree = upgma.run()

        # Printing the tree with a mapping for the species
        print("\nPhylogenetics Tree created:\n")
        tree.print_tree({
            idx: self.get_specie_from_seq(seq)
            for idx, seq in enumerate(best_seqs)
        })
        print()

        print("\n\t:::Step 5 - Creating Graph using UPGMA distance matrix and cut values %s:::\n" % str(self.cut))
        for c in self.cut:
            # Creating Graph from distance matrix with the given cut value
            g = MyGraph.create_from_num_matrix(upgma.dists_mat, c)

            # Printing the Graph and its Metrics
            print("Cut value: %s\n" % str(c))
            g.print_graph_and_metrics()
