from abc import ABC, abstractmethod
import pickle

class Seq(ABC):
    """Abstract class that represents a sequence"""

    class InvalidSequenceException(Exception):
        """Error raised when sequence contains invalid characters"""
        pass

    def __init__(self, seq):
        self._seq = seq.upper()

        if (not self.validate()):
            raise self.InvalidSequenceException()
            

    def get_seq(self):
        return self._seq

    def __len__(self):
        return len(self._seq)

    @abstractmethod
    def validate(self):
        pass

    @abstractmethod
    def __str__(self):
        pass

    def pretty_print(self):
        print("* Sequence:")
        for i in range(0, len(self._seq), 60):
            print(self._seq[i: i + 60])

    def freq(self):
        """Computes the frequency of the symbols in the stored sequence and 
        returns a sorted dictionary"""
        symbols = {}
        for base in self._seq:
            if (base in symbols):
                symbols[base] += 1
            else:
                symbols[base] = 0

        return sorted(symbols.items(),
                      key=lambda symbol: symbol[1],
                      reverse=True)

    @staticmethod
    def readFile(fileName):
        return open(fileName, "r")

    @staticmethod
    def writeFile(fileName):
        return open(fileName, "w+")

    def readSequence(self, fileName):
        self._seq = ""
        for line in Seq.readFile(fileName):
            self._seq += line.strip()

    def writeSequence(self, fileName):
        f = Seq.writeFile(fileName)
        for i in range(0, len(self._seq), 60):
            f.write(self._seq[i: i + 60] + '\n')
        f.close()

    def save(self, fileName):
        bf = open(fileName, mode='wb')
        pickle.dump(self, bf)
        bf.close()
