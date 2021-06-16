# -*- coding: utf-8 -*-

import gzip
import logging
import math
import os
import shutil
from io import IOBase

from src.logging.tqdmLoggingHandler import TqdmLoggingHandler
from tqdm import tqdm


class FastaReader(object):
    """Gets a FASTA file and parse it, getting information about the chromosomes.

    After parse the fasta file, the class will have an attribute (`fasta_data`) with
    information per chromosome in a dictionary data structure type with this keys:

    - **label_length**: length of the chromosome label.
    - **line_start**: the line where the chromosome starts in fasta file.
    - **index_start**: index where the chromosome starts in fasta file.
    - **index_ends**: index where the chromosome ends in fasta file.
    - **index**: index of the chromosome on the local list of chromosomes.
    - **chromosme_length**: number of nucleotides that are in the chromosome.

    For example, after analize a Fasta file with two chromosmes, the dictionary will be:

    ```python
        fasta_data = {
            'chr1': {
                'label_length': ...,
                'line_start': ...,
                'index_start': ...,
                'index_ends': ...,
                'index': ...,
                'chromosme_length': ...,
            },
            'chr2': {
                'label_length': ...,
                'line_start': ...,
                'index_start': ...,
                'index_ends': ...,
                'index': ...,
                'chromosme_length': ...,
            }
        }
    ```

    The class has an attribute (chromosmes) that contains the number and labels of all
    the chrosmomes that are in the fasta file. Using the previous example, this
    attribute would be:

    ```python
        chromosmes = ['chr1', 'chr2']
    ```

    This class allow to get sequences from the fasta file from chromosomes using
    `get_sequence` method.

    Parameters
    ----------
    fasta_path : str
        Path of the fasta file.
    """

    fasta_data: str = {}
    """The dictionary that contains all chromosomes data of the FASTA file, Each key is
    the identifier of each chromosome."""

    chromosmes: list = []
    """List of the identifiers of all the chromosomes of the FASTA file."""

    fasta_filename: str = None
    """Name of the fasta file"""

    _line_length: int = 50
    """Length of the lines that contains nucleotides of the chromosome"""

    def __init__(self, fasta_path: str):
        logging.info("Loading fasta file")
        self.fasta_filename = fasta_path.replace(".gz", "")

        logging.info("Loading fasta information")
        self._set_fasta_file(fasta_path)
        self._get_fasta_data()

    def get_fasta(self) -> IOBase:
        """Returns the fasta file.

        Returns
        -------
        Fasta file.
        """
        return self.fasta_file

    def _set_fasta_line_length(self) -> int:
        """Sets the line length of the fasta file.

        Returns
        -------
        Fasta line length.
        """
        self.fasta_file.seek(0, 0)
        for i in self.fasta_file:
            if not i.startswith(">"):
                self._line_length = len(i) - 1
                break
        return self._line_length

    def _set_fasta_file(self, fasta_path: str):
        """Gets the fasta file from gz file and unzip the file to get the complete fasta
        file. If there is a unziped file with the same name without .gz extension, the
        class use that file instead of unzip the gz one.

        Parameters
        ----------
        fasta_path : str
            Path of the fasta file.
        """
        if not os.path.isfile(self.fasta_filename):
            logging.info("Unzip fasta gz file")
            self.fasta_file = gzip.open(fasta_path, "r")

            with open(self.fasta_filename, "wb") as f_out:
                shutil.copyfileobj(self.fasta_file, f_out)

        self.fasta_file = open(self.fasta_filename, "r")

    def _get_chromosome_index(self, chromosome: str) -> int:
        """Returns the index of the chromosome on the fasta file.

        Parameters
        ----------
        chromosome: str
            Label of the chromosome.

        Returns
        -------
        Chromosome index on fasta file.
        """
        line_start = self.fasta_data[chromosome]["line_start"]

        if line_start <= 0:
            return 0

        chromosome_position = self.chromosmes.index(chromosome)
        previous_chromosmes_labels_length = sum(
            [len(self.chromosmes[i]) + 2 for i in range(chromosome_position)]
        )

        number_new_line_char = (line_start - chromosome_position) * (
            self._line_length + 1
        )

        index = number_new_line_char + previous_chromosmes_labels_length

        # If the last line is not complete the index will be greater that it could be,
        # it's necessary to step back to get the real start
        self.fasta_file.seek(index, 0)
        char = self.fasta_file.read(1)
        while char != ">":
            index -= 1
            self.fasta_file.seek(index, 0)
            char = self.fasta_file.read(1)

        return index

    def _set_chromosme_start_data(self, chromosme: str) -> dict:
        """Sets the data about the start of a chromosme on a fasta file.

        The data of the chromosme is stored as a dictionary in the chromosmes data
        dictionary (`fasta_data`) with the keys:

        - **label_length**: length of the chromosome label.
        - **line_start**: the line where the chromosome starts in fasta file.
        - **index_start**: index where the chromosome starts in fasta file.

        Parameters
        ----------
        chromosme : str
            Chromosome labeled as `line_fasta:>chromosme`, which is obtained from
            `.index` file, generated previously.

        Returns
        -------
        Chromosome start data.
        """
        chromosme_id = chromosme.split(":")
        chromosme_label = chromosme_id[1].rstrip().replace(">", "")

        self.chromosmes.append(chromosme_label)

        self.fasta_data[chromosme_label] = {
            "label_length": len(f">{chromosme_label}\n"),
            "line_start": int(chromosme_id[0]) - 1,
        }

        self.fasta_data[chromosme_label]["index_start"] = self._get_chromosome_index(
            chromosme_label
        )

        return self.fasta_data[chromosme_label]

    def _set_chromosme_end_data(self, chromosme_local_index: int) -> dict:
        """Sets the data about the end of a chromosme in the fasta file. This method has
        to be executed after `_set_chromosme_start_data` to work correctly.

        The data of the chromosme is stored as a dictionary in the chromosmes data
        dictionary (`fasta_data`) with the keys:

        - **index_ends**: index where the chromosome ends in fasta file.
        - **index**: index of the chromosome on the local list of chromosomes.
        - **chromosme_length**: number of nucleotides that are in the chromosome.

        Parameters
        ----------
        chromosme_local_index: int
            Index of the chromosme in the local list.

        Returns
        -------
        Chromosome end data.
        """
        chromosome = self.chromosmes[chromosme_local_index]
        chromosome_data = self.fasta_data[chromosome]

        label_length = self.fasta_data[chromosome]["label_length"]
        index_start = self.fasta_data[chromosome]["index_start"]
        is_last = chromosme_local_index == len(self.chromosmes) - 1

        index_ends = os.stat(self.fasta_filename).st_size
        if not is_last:
            next_chromosome = self.chromosmes[chromosme_local_index + 1]
            index_ends = self.fasta_data[next_chromosome]["index_start"]

        number_of_chars = index_ends - (label_length + index_start)
        # Count the number of '\n (new line character)' present into the chromsome
        newlines_chars = math.ceil(number_of_chars / (self._line_length + 1))
        chromosme_length = number_of_chars - newlines_chars

        chromosome_data["index_ends"] = index_ends
        chromosome_data["index"] = chromosme_local_index
        chromosome_data["chromosme_length"] = chromosme_length

        return chromosome_data

    def _get_fasta_data(self) -> dict:
        """Generates the data of the chromosomes from the fasta file.

        The data of each chromosome is stored in the chromosomes data dictionary
        (`fasta_data`), which have a key per chromsome.

        Each chromosome has as value a dictionary with data about the chromsome in the
        fasta file, with the next keys:

        - **label_length**: length of the chromosome label.
        - **line_start**: the line where the chromosome starts in fasta file.
        - **index_start**: index where the chromosome starts in fasta file.
        - **index_ends**: index where the chromosome ends in fasta file.
        - **index**: index of the chromosome on the local list of chromosomes.
        - **chromosme_length**: number of nucleotides that are in the chromosome.

        For example, after analize a Fasta file with two chromosmes, the dictionary
        will be:

        ```python
            fasta_data = {
                'chr1': {
                    'label_length': ...,
                    'line_start': ...,
                    'index_start': ...,
                    'index_ends': ...,
                    'index': ...,
                    'chromosme_length': ...,
                },
                'chr2': {
                    'label_length': ...,
                    'line_start': ...,
                    'index_start': ...,
                    'index_ends': ...,
                    'index': ...,
                    'chromosme_length': ...,
                }
            }
        ```

        Returns
        -------
        Fasta data.
        """
        logger = logging.getLogger()
        tqdm_out = TqdmLoggingHandler(logger, level=logging.INFO)

        fasta_index_filename = f"{self.fasta_filename}.index"
        self.fasta_data = {}
        self.chromosmes = []
        self._set_fasta_line_length()

        """ TODO: Add compatibility with windows using "type file | findstr /R /C:">" """
        command = (
            f"cat {self.fasta_filename} | grep -n '>' > {fasta_index_filename}"
        )
        os.system(command)
        logging.info("Loading chromosomes")
        with open(fasta_index_filename, "r") as fasta_index_file:
            for i in tqdm(fasta_index_file, file=tqdm_out):
                self._set_chromosme_start_data(i)

        logging.info(f"Loading chromosomes sizes")
        for i in tqdm(range(len(self.chromosmes)), file=tqdm_out):
            self._set_chromosme_end_data(i)

        os.remove(fasta_index_filename)

        return self.fasta_data

    def _get_nucleotide_index(self, pos: int, chromosome: str) -> int:
        """Gets the index of a nucleotide by its position in a chromosome.

        Parameters
        ----------
        pos : int
            Index of the nucletoide on the chromosome.
        chromosome : str
            Label of the chromosome.

        IndexError
            When index is greater tha chromsome length or lower than 0.

        Returns
        -------
        Index of the nucleotide on the fasta file.
        """
        chromosme_length = self.fasta_data[chromosome]["chromosme_length"] - 1

        if pos > chromosme_length or pos < 0:
            raise IndexError(
                f"Invalid index, must be in the interval {0}-{chromosme_length}"
            )

        # It's necessary to taking into account that seek method counts new lines as a
        # character, that's why we add new line characters ('\n')
        num_new_lines = int(pos / self._line_length)
        index_start = self.fasta_data[chromosome]["index_start"]
        label_length = self.fasta_data[chromosome]["label_length"]

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

    def _get_prefix(self, pos: int, length: int, chromosome: str) -> str:
        """Gets the prefix of a given length from a given position in a chromosome.

        Parameters
        ----------
        pos : int
            Position on a chromosome.
        length : int
            Length of the prefix.
        chromosome : str
            Label of the chromosome.

        Returns
        -------
        The prefix.
        """
        starts = pos - length
        if starts <= 0:
            length += starts
            starts = 0

        index = self._get_nucleotide_index(starts, chromosome)

        return self._get_from_interval(index, length)

    def _get_suffix(self, pos: int, length: int, chromosome: str) -> str:
        """Gets the suffix of a given length from a given position in a chromosome.

        Parameters
        ----------
        pos : int
            Position on a chromosome.
        length : int
            Length of the suffix.
        chromosome : str
            Label of the chromosome.

        Returns
        -------
        The suffix.
        """
        pos += 1
        chromosme_length = self.fasta_data[chromosome]["chromosme_length"]
        index = self._get_nucleotide_index(pos, chromosome)

        if pos + length >= chromosme_length:
            length = chromosme_length - pos

        return self._get_from_interval(index, length)

    def _get_nucleotides(self, chromosome: str, pos: int, length: int = 1) -> str:
        """Gets a sequence of nucleotides from a chromosome from a given position on the
        chromosome.

        Parameters
        ----------
        chromosome : str
            The chromosome where the nucleotide is going to be obtained.
        pos : int
            Position of the nucleotide on the chromosome.
        length : int = 1
            Length of the sequence.

        Returns
        -------
        The sequence of nucletoides.
        """
        return self._get_from_interval(
            self._get_nucleotide_index(pos, chromosome), length
        )

    def get_sequence(
        self,
        chromosome: str,
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

        pref = self._get_prefix(pos, from_nuc, chromosome)
        nucleotide = self._get_nucleotides(chromosome, pos, length)
        suff = self._get_suffix(pos + length - 1, to_nuc, chromosome)

        return [pref, nucleotide, suff]
