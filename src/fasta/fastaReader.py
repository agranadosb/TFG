# -*- coding: utf-8 -*-

import gzip
import logging
import math
import os
import shutil

from src.fasta.chromosome import Chromosome
from src.logging.tqdmLoggingHandler import TqdmLoggingHandler
from tqdm import tqdm


class FastaReader(object):
    """Gets a FASTA file and parse it, getting information about the chromosomes.

    After parse the fasta file, the class will have a dictionary with information per
    chromosome.

    For example, after analize a Fasta file with two chromosomes, the dictionary will be:

    ```python
        chromosomes = {
            'chr1': Chromosome(...),
            'chr2': Chromosome(...)
        }
    ```

    This class allow to get sequences from the fasta file from chromosomes using
    `get_sequence` method.

    Parameters
    ----------
    fasta_path : str
        Path of the fasta file.
    """

    chromosomes: str = {}
    """The dictionary that contains all chromosomes data of the FASTA file, Each key is
    the identifier of each chromosome."""

    _chromosomes_list: list = []
    """List of the identifiers of all the chromosomes of the FASTA file."""

    fasta_filename: str = None
    """Name of the fasta file"""

    def __init__(self, fasta_path: str):
        logging.info("Loading fasta file")
        self.fasta_filename = fasta_path.replace(".gz", "")

        logging.info("Loading fasta information")
        self._set_fasta_file(fasta_path)
        self._parse_chromosomes()

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

    def _chromosome_index(
        self, chromosome: str, line_start: int, line_length: int
    ) -> int:
        """Returns the index of the chromosome on the fasta file.

        Parameters
        ----------
        chromosome: str
            Label of the chromosome.
        line_start: int
            Line of the FASTA file where the chromosome starts.
        line_length: int
            Line length of the chromosome.

        Returns
        -------
        Chromosome index on fasta file.
        """
        if line_start <= 0:
            return 0

        chromosome_position = self._chromosomes_list.index(chromosome)
        previous_chromosomes_list_labels_length = sum(
            [len(self._chromosomes_list[i]) + 2 for i in range(chromosome_position)]
        )

        number_new_line_char = (line_start - chromosome_position) * (line_length + 1)

        index = number_new_line_char + previous_chromosomes_list_labels_length

        # If the last line is not complete the index will be greater that it could be,
        # it's necessary to step back to get the real start
        self.fasta_file.seek(index, 0)
        char = self.fasta_file.read(1)
        while char != ">":
            index -= 1
            self.fasta_file.seek(index, 0)
            char = self.fasta_file.read(1)

        return index

    def _line_length(self, chromosome: str) -> int:
        """Returns the line length of a chromsome.
        
        Parameters
        ----------
        chromosome: str
            Name of the chromosome.
        
        Returns
        -------
        Line length.
        """
        def get_char(index):
            self.fasta_file.seek(index, 0)
            return self.fasta_file.read(1)

        def get_index_newline(index):
            char = get_char(index)
            while char != "\n":
                index -= 1
                char = get_char(index)
            return index

        index = os.stat(self.fasta_filename).st_size - 1
        if not len(self._chromosomes_list) == 1:
            chromosome_position = self._chromosomes_list.index(chromosome)
            next_chromosome = self._chromosomes_list[chromosome_position + 1]
            index = self[next_chromosome].index_start - 1

        char = get_char(index)
        if char == "\n":
            index -= 1
            char = get_char(index)

        index = get_index_newline(index) - 1
        new_index = get_index_newline(index)

        return index - new_index

    def _set_chromosome(self, chromosome: str) -> dict:
        """Creates a Chromosome object getting the data from the FASTA file.

        Parameters
        ----------
        chromosome : str
            Chromosome labeled as `line_fasta:>chromosome`, which is obtained from
            `.index` file, generated previously.

        Returns
        -------
        Chromosome start data.
        """
        chromosome_id = chromosome.split(":")
        chromosome_label = chromosome_id[1].rstrip().replace(">", "")

        self._chromosomes_list = [chromosome_label] + self._chromosomes_list

        label_length = len(f">{chromosome_label}\n")
        line_length = self._line_length(chromosome_label)
        index_start = self._chromosome_index(
            chromosome_label, int(chromosome_id[0]) - 1, line_length
        )

        index_ends = os.stat(self.fasta_filename).st_size
        if not len(self._chromosomes_list) == 1:
            index_ends = self.chromosomes[self._chromosomes_list[1]].index_start

        number_of_chars = index_ends - (label_length + index_start)
        # Count the number of '\n (new line character)' present into the chromsome
        newlines_chars = math.ceil(number_of_chars / (line_length + 1))
        length = number_of_chars - newlines_chars

        chromosome = Chromosome(
            self.fasta_file,
            name=chromosome_label,
            line_length=line_length,
            label_length=label_length,
            index_start=index_start,
            length=length,
        )

        self.chromosomes[chromosome_label] = chromosome

        return chromosome

    def _parse_chromosomes(self) -> dict:
        """Create all Cromosmes of he FASTA file and append it to a dictionary with the
        names of chromosomes as keys.

        For example, after analize a Fasta file with two chromosomes, the dictionary
        will be:

        ```python
            chromosomes = {
                'chr1': Chromosome(),
                'chr2': Chromosome()
            }
        ```

        Returns
        -------
        Chromosomes.
        """
        logger = logging.getLogger()
        tqdm_out = TqdmLoggingHandler(logger, level=logging.INFO)

        fasta_index_filename = f"{self.fasta_filename}.index"
        self._chromosomes_list = []
        self.chromosomes = {}

        """ TODO: Add compatibility with windows using "type file | findstr /R /C:">" """
        command = f"cat {self.fasta_filename} | grep -n '>' > {fasta_index_filename}"
        os.system(command)

        logging.info("Loading chromosomes")
        with open(fasta_index_filename, "r") as fasta_index_file:
            lines = fasta_index_file.readlines()

            result = {}
            for i in tqdm(range(len(lines) - 1, -1, -1), file=tqdm_out):
                chromosome = self._set_chromosome(lines[i])

                result[chromosome.name] = {
                    "label_length": chromosome.label_length,
                    "index_start": chromosome.index_start,
                    "length": chromosome.length,
                }

        os.remove(fasta_index_filename)

        return result

    def __getitem__(self, key):
        assert isinstance(key, str)

        return self.chromosomes[key]

    def sequence(
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
        return self[chromosome].sequence(nucleotide, pos, from_nuc, to_nuc)
