# -*- coding: utf-8 -*-

import gzip
import logging
import os
import shutil

from src.logging.tqdmLoggingHandler import TqdmLoggingHandler
from tqdm import tqdm

from vcf import Reader as VcfReader


class VcfMutationsReader(object):
    """Gets information of mutations from vcf and fasta files.

    Using a fasta file and a vcf file gets information and parses data
    from fasta file using the mutations on vcf.

    Parameters
    ----------
    vcf_path : str
        Path of the vcf file
    fasta_path : str
        Path of the fasta file
    """

    def __init__(self, vcf_path: str, fasta_path: str):
        logging.info("Loading vcf file")
        self.vcf_file = VcfReader(open(vcf_path, "r"))

        self.fatsa_data = {}
        self.chromosmes = []
        logging.info("Loading fasta file")
        self.fasta_filename = fasta_path.replace(".gz", "")
        self.fasta_filename_index = fasta_path.replace(".gz", ".index")
        self.fasta_file_line_length = 50

        logging.info("Loading fasta information")
        self.set_fasta_file(fasta_path)
        self.get_fasta_data()

    def get_vcf(self):
        """Returns vcf file

        Returns
        -------
        str
            vcf file
        """
        return self.vcf_file

    def get_fasta(self):
        """Returns fasta file

        Returns
        -------
        str
            fasta file
        """
        return self.fasta_file

    def get_fatsa_data(self):
        """Returns information of chromosomes of fasta file

        Returns
        -------
        dict
            fasta information file
        """
        return self.fatsa_data

    def set_fasta_line_length(self):
        """Sets the line length of the fasta file"""
        for i in self.fasta_file:
            if not i.startswith(">"):
                self.fasta_file_line_length = len(i) - 1
                break

    def set_fasta_file(self, fasta_path):
        """Sets the fsata file

        Parameters
        ----------
        fasta_path : str
            Path of the fasta file
        """
        logging.info("Unzip fasta gz file")
        self.fasta_file = gzip.open(fasta_path, "r")

        with open(self.fasta_filename, "wb") as f_out:
            shutil.copyfileobj(self.fasta_file, f_out)

        self.fasta_file = open(self.fasta_filename, "r")

    def get_chromosome_index(self, chromosome: str):
        """Returns index of the chromosome 'chromosome_label' on the fasta file

        Parameters
        ----------
        chromosome : int
            Label of the chromosome

        Returns
        -------
        str
            chromosome index on fasta file
        """
        line_start = self.fatsa_data[chromosome]["line_start"]

        if line_start <= 0:
            return 0

        chromosome_index = self.chromosmes.index(chromosome)
        chromosme_labels_length = sum(
            [len(self.chromosmes[i]) + 2 for i in range(chromosome_index)]
        )

        new_lines_char = (line_start - chromosome_index) * (
            self.fasta_file_line_length + 1
        )

        chromosome_file_index = new_lines_char + chromosme_labels_length

        return chromosome_file_index

    def set_chromosme_start_data(self, chromosme):
        """Adds chromosome data of 'chromosme'

        Adds chromosome data of 'chromosme':
         - label_length: length of the chromosome label
         - line_start: the line where the chromosome starts in fasta file
         - index_start: index where the chromosome starts in fasta file

        Parameters
        ----------
        chromosme : str
            Chromosome label as '>chromosme:line_fasta'
        """
        chromosme_id = chromosme.split(":")
        chromosme_label = chromosme_id[1].rstrip().replace(">", "")

        self.chromosmes.append(chromosme_label)

        self.fatsa_data[chromosme_label] = {
            "label_length": len(f">{chromosme_label}\n"),
            "line_start": int(chromosme_id[0]) - 1,
        }

        self.fatsa_data[chromosme_label]["index_start"] = self.get_chromosome_index(
            chromosme_label
        )

    def set_chromosme_end_data(self, chromosme_local_index):
        """Adds chromosome information about sizes in fasta file

        Adds chromosome information about sizes in fasta file:
         - index_ends: index where the chromosome ends in fasta file
         - next_label_length: label length of the next chromosome
         - file_ends: shows if the chromosome is the last
         - index: chromosome index
         - chromosme_length: number of nucleotides that are in the chromsome

        Parameters
        ----------
        chromosme_local_index : int
            Index of the chromosme in the local list
        """
        chromosome = self.chromosmes[chromosme_local_index]
        chromosome_information = self.fatsa_data[chromosome]

        if chromosme_local_index < len(self.chromosmes) - 1:
            next_chromosome = self.chromosmes[chromosme_local_index + 1]
            index_ends = self.fatsa_data[next_chromosome]["index_start"]
            next_chromosome_starts = self.fatsa_data[next_chromosome]["line_start"]
            current_chromosome_starts = self.fatsa_data[chromosome]["line_start"]

            chromosme_length = (
                next_chromosome_starts - current_chromosome_starts - 1
            ) * self.fasta_file_line_length

            next_label_length = len(f">{next_chromosome}\n")
            file_ends = False
        else:
            index_ends = os.stat(self.fasta_filename).st_size
            index_start = self.fatsa_data[chromosome]["index_start"]
            label_length = self.fatsa_data[chromosome]["label_length"]
            next_label_length = 0
            file_ends = True

            number_of_elements = (index_ends - label_length) - index_start
            chromosme_length = number_of_elements - int(
                number_of_elements / self.fasta_file_line_length
            )

        chromosome_information["index_ends"] = index_ends
        chromosome_information["next_label_length"] = next_label_length
        chromosome_information["file_ends"] = file_ends
        chromosome_information["index"] = chromosme_local_index
        chromosome_information["chromosme_length"] = chromosme_length

    def get_fasta_data(self):
        """Generates fasta file information about the chromosomes

        Generates fasta file information about the chromosomes:
         - label_length: length of the chromosome label
         - line_start: the line where the chromosome starts in fasta file
         - index_start: index where the chromosome starts in fasta file
         - index_ends: index where the chromosome ends in fasta file
         - next_label_length: label length of the next chromosome
         - file_ends: shows if the chromosome is the last
         - index: chromosome index
        """
        logger = logging.getLogger()
        tqdm_out = TqdmLoggingHandler(logger, level=logging.INFO)

        self.set_fasta_line_length()

        command = (
            f"cat {self.fasta_filename} | grep -n '>' > {self.fasta_filename_index}"
        )
        os.system(command)
        logging.info("Loading chromosomes")
        with open(self.fasta_filename_index, "r") as fasta_index_file:
            for i in tqdm(fasta_index_file, file=tqdm_out):
                self.set_chromosme_start_data(i)

        logging.info(f"Loading chromosomes sizes")
        for i in tqdm(range(len(self.chromosmes)), file=tqdm_out):
            self.set_chromosme_end_data(i)

        os.remove(self.fasta_filename_index)

    def get_nucleotid(self, pos, chromosome):
        """Gets the index of a nucleotide by its position in a chromosome

        Parameters
        ----------
        pos : int
            Index on chromosome 'chromosome'
        chromosome : str
            Label of the chromosome

        IndexError
            When index is greater tha chromsome length or lower
            than 0

        Returns
        -------
        str
            index of the nucleotid in the fasta file
        """
        chromosme_length = self.fatsa_data[chromosome]["chromosme_length"] - 1

        if pos > chromosme_length or pos < 0:
            raise IndexError(
                f"Invalid index, must be in the interval {0}-{chromosme_length}"
            )

        # It's necessary to taking into account that seek method counts new lines as a character
        num_new_lines = int(pos / self.fasta_file_line_length)
        index_start = self.fatsa_data[chromosome]["index_start"]
        label_length = self.fatsa_data[chromosome]["label_length"]

        # Get the position of the nucleotid on the file (1 char is one byte, that's why we use seek)
        return pos + index_start + label_length + num_new_lines

    def get_from_interval(self, from_index: int, length: int):
        """Returns the sequence that starts in 'from_index' and has 'length' length

        Parameters
        ----------
        length : int
            Number of elements from the sequence
        from_index : int
            Start of the sequence

        Returns
        -------
        str
            sequence

        Raises
        ------
        IndexError
            When the sequence from the 'from_index' is over the
            end of the file or is before the start of the file
        """
        self.fasta_file.seek(from_index, 0)
        sequence = ""

        while length != 0:
            searched_sequence = self.fasta_file.read(length).split("\n")
            length = len(searched_sequence) - 1
            sequence += "".join(searched_sequence)

        return sequence.upper()

    def get_prefix(self, pos: int, length: int, chromosome: str):
        """Gets prefix from an index on the fasta file

        Parameters
        ----------
        pos : int
            Index on chromosome 'chromosome'
        length : int
            Length of the prefix
        chromosome : str
            Label of the chromosome

        Returns
        -------
        str
            prefix
        """
        pos_initital = pos - length
        if pos - length <= 0:
            pos_initital = 0
            length += pos - length

        index = self.get_nucleotid(pos_initital, chromosome)

        return self.get_from_interval(index, length)

    def get_suffix(self, pos, length, chromosome):
        """Gets suffix from an index on the fasta file

        Parameters
        ----------
        pos : int
            Index on chromosome 'chromosome'
        length : int
            Length of the suffix
        chromosome : str
            Label of the chromosome

        Returns
        -------
        str
            suffix
        """
        pos += 1
        chromosme_length = self.fatsa_data[chromosome]["chromosme_length"]
        index = self.get_nucleotid(pos, chromosome)

        if pos + length >= chromosme_length:
            length = chromosme_length

        return self.get_from_interval(index, length)

    def get_nucleotide_sequence(self, chromosome: str, pos: int, length: int = 1):
        """Gets the nucleotide (or sequence) from the psoition 'pos' of a chromosme

        Parameters
        ----------
        chromosome : str
            The chromosome where the nucleotide is going to be obtained
        pos : int
            Position of the nucleotide on the chromosome 'chromosome'
        length : str, optional
            If its a sequence of nucleotides what it wants to be obtained,
            indicates the length of that sequence

        Returns
        -------
        str
            nucleotid in the fasta
        """
        return self.get_from_interval(
            self.get_nucleotid(pos, chromosome), length
        )

    def get_sequence(
        self,
        chromosome: str,
        nucleotide: str,
        pos: int,
        from_nuc: int,
        to_nuc: int,
    ):
        """Gets the a fasta's sequence by a position in a chromosome.

        Gets the nucleotide in position 'pos' from the chromosme 'chromosome'
        and, if the parameters 'from_nuc' or 'to_nuc' are given, the prefix and suffix
        with 'from_nuc' and 'to_nuc' length respectively as tuple (prefix, nucleotid, suffix)

        Parameters
        ----------
        chromosome : str
            The chromosome where the nucleotide is going to be obtained
        nucleotide : str
            Label of the nucleotide
        pos : int
            Position of the nucleotide on the chromosome 'chromosome'
        from_nuc : int
            Length of the prefix of the sequence
        to_nuc : int
            Length of the suffix of the sequence

        Returns
        -------
        str
            nucleotid in the fasta
        tuple
            nucleotide with the prefix and the suffix
        """
        length = len(nucleotide)

        pref = self.get_prefix(pos, from_nuc, chromosome)
        nucleotide = self.get_nucleotide_sequence(chromosome, pos, length)
        suff = self.get_suffix(pos + length - 1, to_nuc, chromosome)

        return (pref, nucleotide, suff)
