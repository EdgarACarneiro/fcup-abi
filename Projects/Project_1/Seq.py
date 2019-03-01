from abc import ABC, abstractmethod


class Seq(ABC):
    """Abstract class that represents a sequence"""

    class InvalidSequenceException(Exception):
        """Error raised when sequence contains invalid characters"""
        pass

    def __init__(self, seq):
        self.__seq = seq

    @abstractmethod
    def validate(self):
        pass