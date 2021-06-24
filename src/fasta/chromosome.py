from typing import TextIO


class Chromosome(object):
    """Class that represents a chromosome.

    Using this class a sequence can be obtained as list, for example, if we want to get
    a sequence that starts at position 123 and finish at 456, we can get it using:

    ```python
        chromosome[123:457]
    ```

    If we want to get a prefix from a position, for example, 7 symbols:

    ```python
        chromosome[123:-7]
    ```

    Parameters
    ----------
    fasta_file: TextIO
        Fasta file.
    name: str
        Name of the chromosome.
    line_length: int
        Line length of the chromosome.
    label_length: int
        Length of the first line of the chromosome.
    index_start: int
        Index where the chrosomosme starts on the FASTA file.
    length: int
        Length of the chromosome.
    """

    def __init__(
        self,
        fasta_file: TextIO,
        name: str,
        line_length: int,
        label_length: int,
        index_start: int,
        length: int,
        labels: list = False,
    ) -> None:
        self.name = name
        self.fasta_file = fasta_file
        self.line_length = line_length
        self.label_length = label_length
        self.index_start = index_start
        self.length = length
        self.labels = labels

    def sequence(
        self,
        nucleotide: str,
        pos: int,
        from_nuc: int,
        to_nuc: int,
    ) -> list:
        """Gets a sequence from the fasta file by the position of a nucleotide on a
        chromosome, with a specified prefix length and suffix length.

        Parameters
        ----------
        chromosome : str
            The chromosome where the sequence is going to be obtained.
        nucleotide : str
            Label of the nucleotide.
        pos : int
            Position of the nucleotide on the chromosome.
        from_nuc : int
            Length of the prefix of the sequence.
        to_nuc : int
            Length of the suffix of the sequence.

        Returns
        -------
        The sequence divided in (prefix, nucleotide, suffix).
        """
        length = len(nucleotide)

        pref = self[pos:-from_nuc]
        nucleotide = self[pos : pos + length]
        suff = self[pos + length : pos + length + to_nuc]

        return [pref, nucleotide, suff]

    def _get_nucleotide_index(self, pos: int) -> int:
        """Gets the index of a nucleotide by its position in a chromosome.

        Parameters
        ----------
        pos : int
            Index of the nucletoide on the chromosome.

        IndexError
            When index is greater tha chromsome length or lower than 0.

        Returns
        -------
        Index of the nucleotide on the fasta file.
        """
        length = self.length - 1

        if pos > length or pos < 0:
            raise IndexError(f"Invalid index, must be in the interval {0}-{length}")

        # It's necessary to taking into account that seek method counts new lines as a
        # character, that's why we add new line characters ('\n')
        num_new_lines = int(pos / self.line_length)
        index_start = self.index_start
        label_length = self.label_length

        # Get the position of the nucleotid on the file (1 char is one byte, that's why
        # we use seek)
        return pos + index_start + label_length + num_new_lines

    def _get_from_interval(self, starts: int, length: int) -> str:
        """Returns a sequence that starts at a given index of the fasta file and has a
        given length.

        Parameters
        ----------
        length : int
            Length of the sequence.
        starts : int
            Index of the fasta file where the sequence starts.

        Raises
        ------
        IndexError
            When the start position is wrong or invalid.

        Returns
        -------
        The sequence.
        """
        self.fasta_file.seek(starts, 0)
        sequence = ""

        while length != 0:
            searched_sequence = self.fasta_file.read(length).split("\n")
            # Number of newline characters present on the sequence searched
            length = len(searched_sequence) - 1
            sequence += "".join(searched_sequence)

        return sequence.upper()

    def __getitem__(self, key):
        if isinstance(key, (int)):
            key = slice(key, key + 1, None)

        start = key.start or 0

        if start < 0:
            raise IndexError()

        stop = key.stop or self.length
        if stop < 0:
            stop, start = start, start + stop
            if start < 0:
                start = 0

        if start < 0:
            start = 0
        if stop > self.length:
            stop = self.length

        length = stop - start
        return self._get_from_interval(self._get_nucleotide_index(start), length)

    def __ln__(self):
        return self.length
