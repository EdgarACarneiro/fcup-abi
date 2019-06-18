class SubstMatrix:

    def __init__(self):
        self.alphabet = ""
        self.sm = {}

    def __getitem__(self, ij):
        """Get a cell or row from the substitution matrix"""
        i, j = ij
        return self.score_pair(i, j)

    def score_pair(self, c1, c2):
        """Scores the given pair using the stored substitution matrix"""
        if c1 not in self.alphabet or c2 not in self.alphabet:
            return None
        return self.sm[c1+c2]

    @staticmethod
    def read_submat_file(filename, sep='\t'):
        """Read a substitution matrix from the given file and
        creates a new SubstMatrix object with it."""
        new_sm = SubstMatrix()

        with open(filename, "r") as f:
            new_sm.alphabet = [symbol for symbol
                             in f.readline().replace('\n', '').split(sep)]

            new_sm.sm = {}
            for i, line in enumerate(f):
                line_symbol = line.replace('\n', '').split(sep)

                for j in range(0, len(line_symbol)):
                    new_sm.sm[new_sm.alphabet[i] + new_sm.alphabet[j]] =\
                        int(line_symbol[j])

        return new_sm

    @staticmethod
    def create_submat(alphabet, match, mismatch):
        """Create a new SubstMatrix object with the
        given alphabet using the given values for match
        and mismatch"""
        new_sm = SubstMatrix()

        new_sm.alphabet = alphabet
        new_sm.sm = {i + j: match if i == j else mismatch
                   for i in alphabet for j in alphabet}

        return new_sm
