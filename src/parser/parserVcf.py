# -*- coding: utf-8 -*-

import logging
from abc import ABC, abstractmethod

from src.vcf.vcfReader import VcfMutationsReader


class ParserVcf(ABC):
    """Parses data from vcf and fasta and prepare that data for a
    machine learning model

    Parameters
    ----------
    vcf_path : str
        Path of the vcf file
    fasta_path : str
        Path of the fasta file
    """

    name = False

    def __init__(self, vcf_path: str, fasta_path: str):
        self.vcf_reader = VcfMutationsReader(vcf_path, fasta_path)
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
        VcfMutationsReader
            vcf reader
        """
        return self.vcf_reader

    """ TODO: add retrive_string_sequence abstract method """

    @abstractmethod
    def method(self, sequence: tuple, mutation: str):
        """Generates the result sequence from the original and the mutation

        Parameters
        ----------
        sequence : tuple
            Tuple with the prefix, the nucleotide and the suffix of the sequence
        muatation : str
            Mutation of the nucleotide

        Returns
        -------
        The parsed sequence

        """
        pass

    def original_sequence_to_string(self, prefix, sequence):
        """Gets the original sequence and generates a string representation

        Parameters
        ----------
        prefix : str
            Prefix to append before the original sequence
        sequence:
            Sequence to be representated

        Returns
        -------
        String representation of the original sequence

        """
        return f'{prefix}{"".join(sequence)}\n'

    def default_filename(self):
        """Returns a default filename with the results

        Returns
        -------
        Filename
        """
        return f"parsed_{self.name}_data.pvcf"

    def sequence_to_string(self, original_sequence, prefix, sequence, mutation):
        """Gets a sequence and generates a string representation

        Parameters
        ----------
        original_sequence : str
            String representation of the original sequence
        prefix : str
            Prefix to append before the parsed sequence
        sequence:
            Sequence to be representated
        mutation:
            Mutation of the nucleotide

        Returns
        -------
        String representation of the sequence

        """
        return (
            f'{original_sequence}{prefix}{"".join(self.method(sequence, mutation))}\n'
        )

    def generate_sequences(
        self,
        path: str,
        filename: str = False,
        write_chromosme: bool = False,
        add_original: bool = True,
        prefix_length=5,
        suffix_length=5,
    ):
        """Generates a file with the sequences mutated using 'method'

        Parameters
        ----------
        path : str
            Path to store the data
        filename : str, optional
            Filename of the result file
        write_chromosme : bool, optional
            Indicates the chromosme where the sequences being
        add_original : bool, optional
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
