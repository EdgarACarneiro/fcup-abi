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
        return []

    def add_sequence_database(self, seq):
        """Add an extra sequence to the database"""
        self.database.append(seq)

    def build_map(self, query_seq):
        """Perform the hashing of the query sequence, i.e.
        pre-process the query sequence to identify all words
        of size word_size keeping the positions where they occur"""

        for i in range(len(query_seq) - self.word_size + 1):
            subseq = query_seq[i: i + self.word_size]

            if subseq in self.mapping:
                self.mapping[subseq].append(i)
            else:
                self.mapping[subseq] = [i]

    def get_hits(self, target_seq, threshold=None):
        if threshold is None:
            threshold = self.word_size

        res = []  # list of hits: (index in query, index in target)

        for i in range(len(target_seq) - threshold + 1):
            subseq = target_seq[i: i + threshold]

            # Subseq of target sequence in mapping of query sequence
            if subseq in self.mapping:
                hit = self.mapping[subseq]

                for ind in hit:
                    res.append((ind, i))

        return res

    # TODO Develop a variation of the getHits function that allows
    # for mismatches. It should allow at most 1 mismatch and return
    # all hits that have w or w-1 matches.

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
        hits = self.get_hits(target, query)
        best_align = (None, None, None, -1.0)

        for h in hits:
            align = self.extends_hit(target, h, query)

            # Better score or same score but smaller size
            if align[3] > best_align[3] or \
                    (align[3] == best_align[3] and align[2] < best_align[2]):
                best_align = align

        # Same return tuple as extends_hit
        return best_align if best_align[0] != None else ()

    def best_alignment(self, query):
        self.build_map(query)
        best = (0, 0, 0, -1.0, 0)

        for k in range(0, len(self.database)):
            best_hit_align = self.hit_best_score(self.database[k], query)

            # Better score or same score but smaller size
            if best_hit_align != () and \
                    (best_hit_align[3] > best[3] or
                     (best_hit_align[3] == best[3] and best_hit_align[2] < best[2])):
                best = (*best_hit_align, k)

        return best if best[3] >= 0 else ()

def test1():
    mb = MyBlast("seqBlast.txt", 11)
    query = "gacgcctcgcgctcgcgcgctgaggcaaaaaaaaaaaaaaaaaaaatcggatagctagctgagcgctcgatagcgcgttcgctgcatcgcgtatagcgctgaagctcccggcgagctgtctgtaaatcggatctcatctcgctctatcct"
    r = mb.best_alignment(query)
    print(r)


def test2():
    mb = MyBlast("seqBlast.txt", 11)
    query2 = "cgacgacgacgacgaatgatg"
    r = mb.best_alignment(query2)
    print(r)
