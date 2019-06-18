class MyBlast:
    """Simple BLAST implementation.
    Without substitution matrices.
    Only perfect hits are considered, i.e. the threshold T is equal to w."""

    database: list
    word_size: int
    mapping: dict

    def __init__(self, filename=None, word_size=3):
        self.database =\
            MyBlast.read_database(filename)\
            if filename is not None\
            else []

        self.word_size = word_size
        self.mapping = {}

    @staticmethod
    def read_database(filename):
        """From file with sequences line by line read the sequences to a list"""
        with open(filename, 'r') as fd:
            return [line.strip() for line in fd]
        return []

    @staticmethod
    def read_query(filename):
        """Read a query sequence froma a given file"""
        with open(filename, 'r') as fd:
            for line in fd:
                return line.strip()

    def add_sequence_database(self, seq):
        """Add an extra sequence to the database"""
        self.database.append(seq)

    def build_map(self, query_seq):
        """Perform the hashing of the query sequence, i.e.
        pre-process the query sequence to identify all words
        of size word_size keeping the positions where they occur"""
        self.mapping = {}  # Clearing dictionary

        for i in range(len(query_seq) - self.word_size + 1):
            subseq = query_seq[i: i + self.word_size]

            if subseq in self.mapping:
                self.mapping[subseq].append(i)
            else:
                self.mapping[subseq] = [i]

    def get_hits(self, target_seq):
        res = []  # list of hits: (index in query, index in target)

        for i in range(len(target_seq) - self.word_size + 1):
            subseq = target_seq[i: i + self.word_size]

            # Subseq of target sequence in mapping of query sequence
            if subseq in self.mapping:
                hit = self.mapping[subseq]

                for ind in hit:
                    res.append((ind, i))

        return res

    def extends_hit(self, target, hit, query):
        """Extend the hits found in get_hits(), in both directions while:
        increase in score is >= positions in tht extension / 2"""
        hit_idx_query, hit_idx_target = hit

        # move forward
        k = 0
        fw_score_inc = 0
        best_k = 0

        while 2 * fw_score_inc >= k and \
                hit_idx_query + self.word_size + k < len(query) and \
                hit_idx_target + self.word_size + k < len(target):
            if query[hit_idx_query + self.word_size + k] == \
                    target[hit_idx_target + self.word_size + k]:
                fw_score_inc += 1
                best_k = k + 1
            k += 1
        size = self.word_size + best_k

        # move backwards
        k = 0
        bw_score_inc = 0
        best_k = 0

        while 2 * bw_score_inc >= k and \
                hit_idx_query > k and \
                hit_idx_target > k:
            if query[hit_idx_query - k - 1] == target[hit_idx_target - k - 1]:
                bw_score_inc += 1
                best_k = k+1
            k += 1
        size += best_k

        # Return (idx align on query, idx align on sequence, size align, score)
        return (hit_idx_query - best_k,
                hit_idx_target - best_k,
                size,
                self.word_size + fw_score_inc + bw_score_inc)

    def hit_best_score(self, target, query):
        """Identify the best alignment between the query sequence and the target sequence"""
        hits = self.get_hits(target)
        best_align = (None, None, None, -1.0)

        for h in hits:
            align = self.extends_hit(target, h, query)

            # Better score or same score but smaller size
            if align[3] > best_align[3] or \
                    (align[3] == best_align[3] and align[2] < best_align[2]):
                best_align = align

        # Same return tuple as extends_hit
        return best_align if best_align[0] != None else ()

    def best_alignments(self, query, top=1):
        """Get the best alignments by computing all
        alignments and return only the 'top' ones"""
        self.build_map(query)
        alignments = []

        for k in range(0, len(self.database)):
            best_hit_align = self.hit_best_score(self.database[k], query)

            if best_hit_align != ():
                alignments.append((*best_hit_align, k))

        # Sorting the list first by score and then by size
        # and return 'top' elements
        # Return (idx align on query, idx align on sequence, size align,
        #         score, index of sequence in db)
        return sorted(alignments, reverse=True,
                      key=lambda align: (align[3], align[2]))[:top]
