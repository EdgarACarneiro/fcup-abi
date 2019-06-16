class MyAlign:

    list_seqs: list

    def __init__(self, lseqs, al_type = "protein"):
        self.list_seqs = lseqs
        self.al_type = al_type
    
    def __len__(self): # number of columns
        return len(self.list_seqs[0])
    
    def __getitem__(self, n):
        if type(n) is tuple and len(n) ==2: 
            i, j = n
            return self.list_seqs[i][j]
        elif type(n) is int: return self.list_seqs[n]
        return None
    
    def __str__(self):
        res = ""
        for seq in self.list_seqs:
            res += "\n" + seq 
        return res
    
    def num_seqs(self):
        return len(self.list_seqs)
   
       
    def column (self, index):
        res = []
        for k in range(len(self.list_seqs)):
            res.append(self.list_seqs[k][index])
        return res
    
    def consensus (self):
        return "".join(
            [max(set(
                    filter(lambda el: el != '-', self.column(i))
                ))
                for i
                in range(0, len(self))]
        )


if __name__ == "__main__": 
    alig = MyAlign(["ATGA-A","AA-AT-"], "dna")
    print(alig)
    print(len(alig))
    print(alig.column(2))
    print(alig[1,1])
    print(alig[0])
    print(alig.consensus())

    alig2 = MyAlign(["VJKK","JRSK","VRSK"])
    #print(alig2.richpolarbasic())