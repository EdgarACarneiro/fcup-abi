from Seq import Seq


class Protein(Seq):
    
    def __init__(self, seq):
        super().__init__(seq)

    def validate(self):
        return all(n.isalpha() or n is '_' or n is '*' for n in self._seq)

    #TODO: revise name of functions and compare them to the teacher functions
    def all_proteins_rf():
        """Computes all possible proteins in the stored aminoacid sequence
        Complexity: O(log(n))"""
        proteins = []
        curr = []
        begin = False
        
        for aa in self._seq:
            if aa is "M":
                curr.append("M")
                begin = True
                continue
            
            if begin:
                curr[len(curr) - 1] += aa
                
            if aa is "_" and begin:
                for i in range(0, len(curr)):
                    seq = ""
                    for j in range(i, len(curr)) :
                        seq += curr[j]
                    proteins.append(seq)
                curr = []
                begin = False
        
        return proteins

def main():
    assert Protein("aomademaedoaed_*adiaedae").validate(), "Invalid protein"

if __name__ == "__main__":
    main()