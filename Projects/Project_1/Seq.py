from abc import ABC, abstractmethod


class Seq(ABC):
    """Abstract class that represents a sequence"""

    class InvalidSequenceException(Exception):
        """Error raised when sequence contains invalid characters"""
        pass

    def __init__(self, seq):
        self._seq = seq.upper()

    def get_seq(self):
        return bio._seq

    @abstractmethod
    def validate(self):
        pass

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

    def readFile(fileName):
        return open(fileName, "r")

    def writeFile(fileName):
        return open(fileName, "w+")

    def readSequence(fileName):
        self._seq = ""
        for line in readFile(fileName):
            self._seq += line.strip()

    def writeSequence(fileName):
        f = writeFile(fileName)
        for i in range(0, len(self._seq), 60):
            f.write(self._seq[i: i + 60] + '\n')
        f.close()