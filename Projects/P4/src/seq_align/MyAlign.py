from seq_align import recover_global_align_multiple_solutions,\
    global_align_multiple_solutions


class MyAlign:

    list_seqs: list
    align_type: str

    def __init__(self, lseqs, al_type='str'):
        """Verifies if all the sequences are from the same type,
        and if so saves them"""
        seq_type = None

        for seq in lseqs:
            if not seq_type:
                seq_type = type(seq)

            if type(seq) != seq_type:
                raise Exception(
                    "All sequences in list must be of the same type")

        if seq_type == str:
            self.list_seqs = lseqs
            self.align_type = al_type
        else:
            self.list_seqs = [str(seq) for seq in lseqs]
            self.align_type = seq_type.__name__

    @staticmethod
    def align_from_global_alignment(seq1, seq2, sm, g):
        """Create MyAlign object from global alignment"""
        return MyAlign(
            recover_global_align_multiple_solutions(
                global_align_multiple_solutions(
                    str(seq1),
                    str(seq2),
                    sm,
                    g)[1],
                str(seq1),
                str(seq2)
            )[0], type(seq1).__name__)

    def __len__(self):
        """Returns the number of columns"""
        return len(self.list_seqs[0])

    def __getitem__(self, n):
        """Gets the cell of the alignment matrix if a tuple is given
        or a row in the matrix (sequence) if a int is given"""
        if type(n) is tuple and len(n) == 2:
            i, j = n
            return self.list_seqs[i][j]

        elif type(n) is int:
            return self.list_seqs[n]

        return None

    def __str__(self):
        """Returns list of sequences separated by '\n'"""
        return "".join(["\n%s" % seq for seq in self.list_seqs])

    def pretty_print(self):
        """Pretty print of the Alignment"""
        for idx, seq in enumerate(self.list_seqs):
            print('%d. %s\n' % (idx, seq))

    def get_seqs(self):
        """Get the list of sequences"""
        return self.list_seqs

    def get_align_type(self):
        """Get the alignment type (DNA, RNA or Protein expected)"""
        return self.align_type

    def num_seqs(self):
        """Returns the number of sequences used in the alignment"""
        return len(self.list_seqs)

    def column(self, index):
        """Returns list with symbols in column index"""
        return [self.list_seqs[k][index]
                for k in range(len(self.list_seqs))]

    def consensus(self):
        """Gets the consensus sequence.
        A consensus sequence is the sequence composed of the most frequent
        character in each column of the alignment ignoring the gap"""
        return "".join(
            [max(set(
                filter(lambda el: el != '-', self.column(i))
            ))
                for i
                in range(0, len(self))]
        )
