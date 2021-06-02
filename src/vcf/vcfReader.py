# -*- coding: utf-8 -*-

import gzip
import logging
import math
import os
import shutil
from io import IOBase

from src.logging.tqdmLoggingHandler import TqdmLoggingHandler
from tqdm import tqdm

from vcf import Reader as VcfReader


class FastaReader(object):
    """Gets vcf and fasta file and generates data about the fasta file and the
    chromosomes that are in the fasta file.

    After parse the fasta file data, the class will have an attribute (fasta_data) with
    information per chromosome in a dictionary data structure:

    - **label_length**: length of the chromosome label
    - **line_start**: the line where the chromosome starts in fasta file
    - **index_start**: index where the chromosome starts in fasta file
    - **index_ends**: index where the chromosome ends in fasta file
    - **index**: chromosome index on the local list of chromosomes
    - **chromosme_length**: number of nucleotides that are in the chromosome

    For example, after analize a Fasta file with two chromosmes, the dictionary will be:

        {
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

    The class has an attribute (chromosmes) that contains the number and labels of all
    the chrosmomes that are in the fasta file. Getting the previous example, this
    attribute would be:

        chromosmes = ['chr1', 'chr2']

    This class allow to get sequences from the fasta file and from chrosomes using the
    methods:
         - get_sequence
         - get_prefix
         - get_nucleotides
         - get_suffix

    Parameters
    ----------
    vcf_path : str
        Path of the vcf file
    fasta_path : str
        Path of the fasta file
    """

    vcf_file: VcfReader = None
    fasta_data: str = {}
    chromosmes: list = []
    fasta_filename: str = None
    fasta_index_filename: str = None
    line_length: int = 50

    def __init__(self, vcf_path: str, fasta_path: str):
        logging.info("Loading vcf file")
        self.vcf_file = VcfReader(open(vcf_path, "r"))

        logging.info("Loading fasta file")
        self.fasta_filename = fasta_path.replace(".gz", "")
        self.fasta_index_filename = fasta_path.replace(".gz", ".index")

        logging.info("Loading fasta information")
        self.set_fasta_file(fasta_path)
        self.get_fasta_data()

    def get_vcf(self) -> VcfReader:
        """Returns vcf file

        Returns
        -------
        Vcf file
        """
        return self.vcf_file

    def get_fasta(self) -> IOBase:
        """Returns fasta file

        Returns
        -------
        Fasta file
        """
        return self.fasta_file

    def set_fasta_line_length(self) -> int:
        """Sets the line length of the fasta file

        Returns
        -------
        Fasta line length
        """
        self.fasta_file.seek(0, 0)
        for i in self.fasta_file:
            if not i.startswith(">"):
                self.line_length = len(i) - 1
                break
        return self.line_length

    def set_fasta_file(self, fasta_path: str):
        """Gets the fasta file from gz file and unzip the file to get the complete fasta
        file.

        Parameters
        ----------
        fasta_path : str
            Path of the fasta file
        """
        if not os.path.isfile(self.fasta_filename):
            logging.info("Unzip fasta gz file")
            self.fasta_file = gzip.open(fasta_path, "r")

            with open(self.fasta_filename, "wb") as f_out:
                shutil.copyfileobj(self.fasta_file, f_out)

        self.fasta_file = open(self.fasta_filename, "r")

    def get_chromosome_index(self, chromosome: str) -> int:
        """Returns the index of the chromosome on the fasta file

        Parameters
        ----------
        chromosome: str
            Label of the chromosome

        Returns
        -------
        Chromosome index on fasta file
        """
        line_start = self.fasta_data[chromosome]["line_start"]

        if line_start <= 0:
            return 0

        chromosome_position = self.chromosmes.index(chromosome)
        previous_chromosmes_labels_length = sum(
            [len(self.chromosmes[i]) + 2 for i in range(chromosome_position)]
        )

        number_new_line_char = (line_start - chromosome_position) * (
            self.line_length + 1
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

    def set_chromosme_start_data(self, chromosme: str) -> dict:
        """Set the data about the start of a chromosme on a fasta file

        The data of the chromosme is stored as a dictionary in the chromosmes data
        dictionary with the keys:

         - **label_length**: length of the chromosome label
         - **line_start**: the line where the chromosome starts in fasta file
         - **index_start**: index where the chromosome starts in fasta file

        Parameters
        ----------
        chromosme : str
            Chromosome label as 'line_fasta:>chromosme', which is obtained from .index
            generated file

        Returns
        -------
        Chromosome start data
        """
        chromosme_id = chromosme.split(":")
        chromosme_label = chromosme_id[1].rstrip().replace(">", "")

        self.chromosmes.append(chromosme_label)

        self.fasta_data[chromosme_label] = {
            "label_length": len(f">{chromosme_label}\n"),
            "line_start": int(chromosme_id[0]) - 1,
        }

        self.fasta_data[chromosme_label]["index_start"] = self.get_chromosome_index(
            chromosme_label
        )

        return self.fasta_data[chromosme_label]

    def set_chromosme_end_data(self, chromosme_local_index: int) -> dict:
        """Set the data about the end of a chromosme in a fasta file. This method has to
        be executed after set_chromosme_start_data to work correctly.

        The data of the chromosme is stored as a dictionary in the chromosmes data
        dictionary with the keys:

         - **index_ends**: index where the chromosome ends in fasta file
         - **index**: chromosome index on the local list of chromosomes
         - **chromosme_length**: number of nucleotides that are in the chromosome

        Parameters
        ----------
        chromosme_local_index: int
            Index of the chromosme in the local list

        Returns
        -------
        Chromosome end data
        """
        chromosome = self.chromosmes[chromosme_local_index]
        chromosome_data = self.fasta_data[chromosome]

        label_length = self.fasta_data[chromosome]["label_length"]
        index_start = self.fasta_data[chromosome]["index_start"]

        if not chromosme_local_index == len(self.chromosmes) - 1:
            next_chromosome = self.chromosmes[chromosme_local_index + 1]
            index_ends = self.fasta_data[next_chromosome]["index_start"]
        else:
            index_ends = os.stat(self.fasta_filename).st_size

        number_of_chars = index_ends - (label_length + index_start)
        newlines_chars = math.ceil(number_of_chars / (self.line_length + 1))
        chromosme_length = number_of_chars - newlines_chars

        chromosome_data["index_ends"] = index_ends
        chromosome_data["index"] = chromosme_local_index
        chromosome_data["chromosme_length"] = chromosme_length

        return chromosome_data

    def get_fasta_data(self) -> dict:
        """Generates data about the chromosomes in the fasta file.

        The data of each chromosme is stored as a dictionary in the chromosmes data
        dictionary with the keys:

        - **label_length**: length of the chromosome label
        - **line_start**: the line where the chromosome starts in fasta file
        - **index_start**: index where the chromosome starts in fasta file
        - **index_ends**: index where the chromosome ends in fasta file
        - **index**: chromosome index on the local list of chromosomes
        - **chromosme_length**: number of nucleotides that are in the chromosome

        Returns
        -------
        Fasta data
        """
        logger = logging.getLogger()
        tqdm_out = TqdmLoggingHandler(logger, level=logging.INFO)

        self.fasta_data = {}
        self.chromosmes = []
        self.set_fasta_line_length()

        """ TODO: Add compatibility with windows using "type file | findstr /R /C:">" """
        command = (
            f"cat {self.fasta_filename} | grep -n '>' > {self.fasta_index_filename}"
        )
        os.system(command)
        logging.info("Loading chromosomes")
        with open(self.fasta_index_filename, "r") as fasta_index_file:
            for i in tqdm(fasta_index_file, file=tqdm_out):
                self.set_chromosme_start_data(i)

        logging.info(f"Loading chromosomes sizes")
        for i in tqdm(range(len(self.chromosmes)), file=tqdm_out):
            self.set_chromosme_end_data(i)

        os.remove(self.fasta_index_filename)

        return self.fasta_data

    def get_nucleotide_index(self, pos: int, chromosome: str) -> int:
        """Gets the index of a nucleotide by its position in a chromosome

        Parameters
        ----------
        pos : int
            Index of the nucletoide on the chromoome
        chromosome : str
            Label of the chromosome

        IndexError
            When index is greater tha chromsome length or lower than 0

        Returns
        -------
        Index of the nucleotide on the fasta file
        """
        chromosme_length = self.fasta_data[chromosome]["chromosme_length"] - 1

        if pos > chromosme_length or pos < 0:
            raise IndexError(
                f"Invalid index, must be in the interval {0}-{chromosme_length}"
            )

        # It's necessary to taking into account that seek method counts new lines as a
        # character
        num_new_lines = int(pos / self.line_length)
        index_start = self.fasta_data[chromosome]["index_start"]
        label_length = self.fasta_data[chromosome]["label_length"]

        # Get the position of the nucleotid on the file (1 char is one byte, that's why
        # we use seek)
        return pos + index_start + label_length + num_new_lines

    def get_from_interval(self, starts: int, length: int) -> str:
        """Returns a sequence that starts at a given index of the fasta file and has a
        given length.

        Parameters
        ----------
        length : int
            Length of the sequence
        starts : int
            Index of the fasta file where the sequence starts

        Raises
        ------
        IndexError
            When the start position is wrong or invalid

        Returns
        -------
        The sequence
        """
        self.fasta_file.seek(starts, 0)
        sequence = ""

        while length != 0:
            searched_sequence = self.fasta_file.read(length).split("\n")
            # Number of newline characters present on the sequence searched
            length = len(searched_sequence) - 1
            sequence += "".join(searched_sequence)

        return sequence.upper()

    def get_prefix(self, pos: int, length: int, chromosome: str) -> str:
        """Gets the prefix of a given length from a given position in a chromosome

        Parameters
        ----------
        pos : int
            Position on a chromosome
        length : int
            Length of the prefix
        chromosome : str
            Label of the chromosome

        Returns
        -------
        The prefix
        """
        starts = pos - length
        if starts <= 0:
            length += starts
            starts = 0

        index = self.get_nucleotide_index(starts, chromosome)

        return self.get_from_interval(index, length)

    def get_suffix(self, pos: int, length: int, chromosome: str) -> str:
        """Gets the prefix of a given length from a given position in a chromosome

        Parameters
        ----------
        pos : int
            Position on a chromosome
        length : int
            Length of the suffix
        chromosome : str
            Label of the chromosome

        Returns
        -------
        The suffix
        """
        pos += 1
        chromosme_length = self.fasta_data[chromosome]["chromosme_length"]
        index = self.get_nucleotide_index(pos, chromosome)

        if pos + length >= chromosme_length:
            length = chromosme_length - pos

        return self.get_from_interval(index, length)

    def get_nucleotides(self, chromosome: str, pos: int, length: int = 1) -> str:
        """Gets a sequence of nucletoides from a chrosome from a given position on the
        chromosome. The default length of the sequence is 1.

        Parameters
        ----------
        chromosome : str
            The chromosome where the nucleotide is going to be obtained
        pos : int
            Position of the nucleotide on the chromosome
        length : int = 1
            Length of the sequence

        Returns
        -------
        The sequence of nucletoides
        """
        return self.get_from_interval(
            self.get_nucleotide_index(pos, chromosome), length
        )

    def get_sequence(
        self,
        chromosome: str,
        nucleotide: str,
        pos: int,
        from_nuc: int,
        to_nuc: int,
    ) -> list:
        """Gets a fasta's sequence by the position on a chromosome, with a given prefix
        and a suffix length.

        Parameters
        ----------
        chromosome : str
            The chromosome where the sequence is going to be obtained
        nucleotide : str
            Label of the nucleotide
        pos : int
            Position of the nucleotide on the chromosome
        from_nuc : int
            Length of the prefix of the sequence
        to_nuc : int
            Length of the suffix of the sequence

        Returns
        -------
        The sequence divided in (prefix, nucleotide, suffix)
        """
        length = len(nucleotide)

        pref = self.get_prefix(pos, from_nuc, chromosome)
        nucleotide = self.get_nucleotides(chromosome, pos, length)
        suff = self.get_suffix(pos + length - 1, to_nuc, chromosome)

        return [pref, nucleotide, suff]
