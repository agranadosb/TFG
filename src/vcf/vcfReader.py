# -*- coding: utf-8 -*-

import gzip
import os
import shutil

from vcf import Reader as VcfReader


class VcfMutationsReader(object):
    """Gets information of mutations by a vcf and a fasta files

    Using a fasta file and a vcf file gets information and parses data
    from fasta file using the mutations on vcf

    Parameters
    ----------
    vcf_path : str
        Path of the vcf file
    fasta_path : str
        Path of the fasta file
    """

    def __init__(self, vcf_path: str, fasta_path: str):
        self.vcf_file = VcfReader(open(vcf_path, "r"))

        self.fatsa_chromosme_information = {}
        self.fasta_keys = []
        self.fasta_filename = fasta_path.replace(".gz", "")
        self.fasta_filename_index = fasta_path.replace(".gz", ".index")
        self.fasta_file_line_length = 50

        self.set_fasta_file(fasta_path)
        self.generate_fasta_information()

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

    def get_fatsa_chromosme_information(self):
        """Returns information of chromosomes of fasta file

        Returns
        -------
        dict
            fasta information file
        """
        return self.fatsa_chromosme_information

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
        self.fasta_file = gzip.open(fasta_path, "r")

        with open(self.fasta_filename, "wb") as f_out:
            shutil.copyfileobj(self.fasta_file, f_out)

        self.fasta_file = open(self.fasta_filename, "r")

    def set_chromosme_information(self, chromosme):
        """Adds chromosome information of 'chromosme'

        Adds chromosome information of 'chromosme':
         - label_length: length of the chromosome label
         - line_starts: the line where the chromosome starts in fasta file
         - index_starts: index where the chromosome starts in fasta file

        Parameters
        ----------
        chromosme : str
            Chromosome label as '>chromosme:line_fasta'
        """
        index_id = chromosme.split(":")
        chromosme_label = index_id[1].rstrip().replace(">", "")

        self.fasta_keys.append(chromosme_label)

        self.fatsa_chromosme_information[chromosme_label] = {
            "label_length": len(f">{chromosme_label}"),
            "line_starts": int(index_id[0]) - 1,
        }

        self.fatsa_chromosme_information[chromosme_label][
            "index_starts"
        ] = self.get_chromosome_file_index(chromosme_label)

    def set_chromosme_sizes(self, chromosme_local_index):
        """Adds chromosome information about sizes in fasta file

        Adds chromosome information about sizes in fasta file:
         - index_ends: index where the chromosome ends in fasta file
         - next_chromosome_length: label length of the next chromosome
         - file_ends: shows if the chromosome is the last
         - index: chromosome index

        Parameters
        ----------
        chromosme_local_index : int
            Index of the chromosme in the local list
        """
        chromosome_information = self.fatsa_chromosme_information[
            self.fasta_keys[chromosme_local_index]
        ]
        index_ends = os.stat(self.fasta_filename).st_size
        next_chromosome_length = 0
        file_ends = True

        if chromosme_local_index != len(self.fasta_keys) - 1:
            next_chromosome = self.fasta_keys[chromosme_local_index + 1]
            index_ends = self.fatsa_chromosme_information[next_chromosome][
                "index_starts"
            ]
            next_chromosome_length = len(f">{next_chromosome}")
            file_ends = False

        chromosome_information["index_ends"] = index_ends
        chromosome_information["next_chromosome_length"] = next_chromosome_length
        chromosome_information["file_ends"] = file_ends
        chromosome_information["index"] = chromosme_local_index

    """ TODO: Añadir el fasta information a una clase a parte """

    def generate_fasta_information(self):
        """Generates fasta file information about the chromosomes

        Generates fasta file information about the chromosomes:
         - label_length: length of the chromosome label
         - line_starts: the line where the chromosome starts in fasta file
         - index_starts: index where the chromosome starts in fasta file
         - index_ends: index where the chromosome ends in fasta file
         - next_chromosome_length: label length of the next chromosome
         - file_ends: shows if the chromosome is the last
         - index: chromosome index
        """
        self.set_fasta_line_length()

        command = (
            f"cat {self.fasta_filename} | grep -n '>' > {self.fasta_filename_index}"
        )
        os.system(command)
        with open(self.fasta_filename_index, "r") as fasta_index_file:
            for i in fasta_index_file:
                self.set_chromosme_information(i)

        for i in range(len(self.fasta_keys)):
            self.set_chromosme_sizes(i)

    def get_nucleotid_fasta_index(self, pos, chromosome):
        """Gets the index of a nucleotide by its position in a chromosome

        Parameters
        ----------
        pos : int
            Index on chromosome 'chromosome'
        chromosome : str
            Label of the chromosome

        Returns
        -------
        str
            index of the nucleotid in the fasta file
        """
        # It's necessary to taking into account that seek method counts new lines as a character
        num_new_lines = int(pos / self.fasta_file_line_length)
        index_starts = self.fatsa_chromosme_information[chromosome]["index_starts"]
        label_length = self.fatsa_chromosme_information[chromosome]["label_length"] + 1

        # Get the position of the nucleotid on the file (1 char is one byte, that's why we use seek)
        return pos + index_starts + label_length + num_new_lines

    def get_sequence_interval(self, from_index: int, length: int):
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

    def get_chromosome_file_index(self, chromosome: str):
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
        line_starts = self.fatsa_chromosme_information[chromosome]["line_starts"]

        if line_starts <= 0:
            return 0

        chromosome_index = self.fasta_keys.index(chromosome)
        chromosme_labels_length = sum(
            [len(self.fasta_keys[i]) + 2 for i in range(chromosome_index)]
        )

        new_lines_char = (line_starts - chromosome_index) * (
            self.fasta_file_line_length + 1
        )

        chromosome_file_index = new_lines_char + chromosme_labels_length

        return chromosome_file_index

    def get_sequence_prefix(self, pos: int, length: int, chromosome: str):
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
        index_starts = self.fatsa_chromosme_information[chromosome]["index_starts"]
        chromosome_label_length = self.fatsa_chromosme_information[chromosome][
            "label_length"
        ]
        index_ends = self.fatsa_chromosme_information[chromosome]["index_ends"]

        last_line = pos % self.fasta_file_line_length == 0
        pointer_index = self.get_nucleotid_fasta_index(pos, chromosome)
        chromosome_starts = index_starts + chromosome_label_length
        sequence_ends = pointer_index - length - 1 + last_line
        interval = sequence_ends - chromosome_starts

        if index_ends - pointer_index <= self.fasta_file_line_length and not last_line:
            sequence_ends += 1

        if interval <= chromosome_label_length and interval <= 0:
            length = interval + chromosome_label_length
            sequence_ends = chromosome_starts + 1

        return self.get_sequence_interval(sequence_ends, length)

    def get_sequence_sufix(self, pos, length, chromosome):
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
        file_ends = self.fatsa_chromosme_information[chromosome]["file_ends"]
        next_chromosome_length = self.fatsa_chromosme_information[chromosome][
            "next_chromosome_length"
        ]
        index_ends = self.fatsa_chromosme_information[chromosome]["index_ends"]

        pointer_index = self.get_nucleotid_fasta_index(pos, chromosome)
        ends = pointer_index + 1 + length

        if ends >= index_ends:
            length_copy = length
            length = index_ends - ends + next_chromosome_length - 1
            if length <= 0 and file_ends:
                length += length_copy

        return self.get_sequence_interval(pointer_index + 1, length)

    def get_nucleotide(self, chromosome: str, pos: int, length: int = 1):
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
        return self.get_sequence_interval(
            self.get_nucleotid_fasta_index(pos, chromosome), length
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

        pref = self.get_sequence_prefix(pos, from_nuc, chromosome)
        nucleotide = self.get_nucleotide(chromosome, pos, length)
        suff = self.get_sequence_sufix(pos + length - 1, to_nuc, chromosome)

        return (pref, nucleotide, suff)
