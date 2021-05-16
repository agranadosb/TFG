# -*- coding: utf-8 -*-

import logging
from abc import ABC, abstractmethod

from src.vcf.vcfReader import FastaReader


class ParserVcf(ABC):
    """Parses data from vcf and fasta and prepare that data for a machine learning model

    Parameters
    ----------
    vcf_path : str
        Path of the vcf file
    fasta_path : str
        Path of the fasta file
    """

    name = False

    def __init__(self, vcf_path: str, fasta_path: str):
        self.vcf_reader = FastaReader(vcf_path, fasta_path)
        logging.info("Loading finalized\n")

    def get_vcf(self):
        """Returns vcf file

        Returns
        -------
        str
            vcf file
        """
        return self.vcf_reader.get_vcf()

    def get_vcf_reader(self):
        """Returns vcf reader

        Returns
        -------
        FastaReader
            vcf reader
        """
        return self.vcf_reader

    @staticmethod
    @abstractmethod
    def retrive_sequence(sequence):
        """Gets a string sequence and returns the sequence in a tuple type

        Parameters
        ----------
        sequence: str
            Sequence in string format

        Returns
        -------
        Sequence in a list format
        """
        pass

    @staticmethod
    @abstractmethod
    def retrive_string_sequence(sequence):
        """Gets a string sequence and returns the sequence in a tuple type

        Parameters
        ----------
        sequence: str
            Sequence in string format

        Returns
        -------
        Sequence in a string format
        """
        pass

    @abstractmethod
    def method(self, sequence: tuple, mutation: str):
        """Parse the sequence with a given mutation to a new sequence with different
        noation, for example:
            from ("a", "b", "c") with mutation = "dd" to ("s", "mm", "s")

        Parameters
        ----------
        sequence : tuple
            Sequence in a tuple shape
        muatation : str
            Mutation of the nucleotide

        Returns
        -------
        The parsed sequence
        """
        pass

    def original_sequence_to_string(self, prefix, sequence):
        """Generates a string of a sequence by append the prefix, infix and suffix and
        append at the starts a given prefix.

        Parameters
        ----------
        prefix : str
            Prefix to append before the original sequence
        sequence:
            Sequence

        Returns
        -------
        Sequence in a string shape with a given prefix

        """
        return f'{prefix}{"".join(sequence)}\n'

    def default_filename(self):
        """Returns a default filename for the results

        Returns
        -------
        Results default filename
        """
        return f"parsed_{self.name}_data.pvcf"

    @abstractmethod
    def sequence_to_string(self, original_sequence, prefix, sequence, mutation):
        """Creates a string from a sequence and a mutation with the original sequence
        (in a string shape) and a given prefix.

        Parameters
        ----------
        original_sequence : str
            Original sequence in a string shape
        prefix : str
            Prefix to append before the parsed sequence
        sequence:
            Sequence
        mutation:
            Mutation

        Returns
        -------
        Sequence and original sequence with the prefix
        """
        pass

    def generate_sequences(
        self,
        path: str,
        filename: str = False,
        write_chromosme: bool = False,
        add_original: bool = True,
        prefix_length=5,
        suffix_length=5,
    ):
        """Generates a file with the sequences mutated using a defined method for parse
        the sequence and the mutation

        Parameters
        ----------
        path : str
            Path to store the data
        [filename : str = default_filename]
            Filename of the result file
        [write_chromosme : bool = False]
            If true adds the chromosme where the sequences being into the file
        [add_original : bool = False]
            If true adds the original sequence into the file
        """
        if not filename:
            filename = self.default_filename()

        with open(f"{path}/{filename}", "w") as parsed_data_file:
            for i in self.get_vcf():
                sequence = self.vcf_reader.get_sequence(
                    i.CHROM, i.REF, i.POS, prefix_length, suffix_length
                )

                prefix = ""
                if write_chromosme:
                    prefix = f"{i.CHROM}\t"

                original_sequence = ""
                if add_original:
                    original_sequence = self.original_sequence_to_string(
                        prefix, sequence
                    )

                parsed_data_file.write(
                    self.sequence_to_string(original_sequence, prefix, sequence, i.ALT)
                )
