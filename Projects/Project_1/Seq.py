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
        """Validates the sequence used in the object creation"""
        pass

    @abstractmethod
    def __str__(self):
        pass

    def pretty_print(self):
        """Method to be overriden by child classes.
        Pretty prints the stored bio sequence"""
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
        """Returns a file descriptor in the read mode for the given file"""
        return open(fileName, "r")

    @staticmethod
    def writeFile(fileName):
        """Returns a file descriptor in the write mode for the given file"""
        return open(fileName, "w+")

    def readSequence(self, fileName):
        """Reads a bio sequence from the given file"""
        self._seq = ""
        for line in Seq.readFile(fileName):
            self._seq += line.strip()

    def writeSequence(self, fileName):
        """Writes the currently stored bio sequence to the given file"""
        f = Seq.writeFile(fileName)
        for i in range(0, len(self._seq), 60):
            f.write(self._seq[i: i + 60] + '\n')
        f.close()

    def save(self, fileName):
        """Saves the bio sequence instance into the given file"""
        bf = open(fileName, mode='wb')
        pickle.dump(self, bf)
        bf.close()
