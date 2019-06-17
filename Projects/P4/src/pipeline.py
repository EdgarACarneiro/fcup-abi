from bioseq import BioSeq, Seq, Dna, Rna, Protein
from seq_align import MyBlast
import re


class Pipeline:

    query_seq: Seq
    query_id: str
    database: {str: Seq}

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

        self.database = db_fasta

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
        return (re.search(r'.*?\[(.*?)\].*', seq_id)).group(1)

    def get_specie(self):
        """Get the query sequence's specie"""
        return Pipeline.get_specie_from_id(self.query_id)

    def execute(self):
        """Execute the Pipeline"""
        # Create copy database without similar species
        db_copy = {_id: seq
                   for _id, seq in self.database.items()
                   if Pipeline.get_specie_from_id(_id) !=
                   self.get_specie()}
        
        # Creating Blast and populating its database
        blast = MyBlast()
        for seq in db_copy.values():
            blast.add_sequence_database(str(seq))

        # Top 10 Alignments
        top_alignments = blast.best_alignments(str(self.query_seq), 10)
        print(top_alignments)




if __name__ == '__main__':
    p = Pipeline('tests/files/source.fasta', 'tests/files/seqdump.txt')
    p.execute()
    # print(p.database)
