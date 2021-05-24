# -*- coding: utf-8 -*-

import logging
from abc import ABC, abstractmethod
from typing import Union

from src.argumentParser.abstractArguments import AbstractParserArguments
from src.vcf.vcfReader import FastaReader
from vcf import Reader as VcfReader


class ParserVcf(AbstractParserArguments, ABC):
    """Parses data from vcf and fasta files and prepare that data for a machine learning
    model.

    The parser needs to define a method named 'method', which generates a new sequence
    from the original one.

    The parser will have static methods to send and retrieve the parsed sequence into a
    file in a specified format. To do this, the parser will need to be defined the
    metohds:

     - **sequence_to_string**: Transforms a sequence in a list or tuple shape into a
    string shape to be written into a file
     - **retrive_sequence**: Get a sequence in a string shape and return as a tuple
     - **retrive_string_sequence**: Get a sequencein a string shape and return as a
    string

    Parameters
    ----------
    vcf_path : str
        Path of the vcf file
    fasta_path : str
        Path of the fasta file
    """

    arguments: list = [
        {
            "key": "amto",
            "name": "add-mutation-to-original",
            "help": "Add mutation to original sequence on parser file",
            "function_argumemnt": {"amto": "add_mutation_to_original"},
            "action": "store_true",
        },
        {
            "key": "ao",
            "name": "add-original",
            "help": "Writes original sequence into parsed sequences file on parser file",
            "function_argumemnt": {"ao": "add_original"},
            "action": "store_true",
        },
        {
            "key": "wc",
            "name": "write-chromosme",
            "help": "Writes the chromsome where the sequence are from on parser file",
            "function_argumemnt": {"wc": "write_chromosme"},
            "action": "store_true",
        },
        {
            "key": "pfilename",
            "name": "parser-filename",
            "help": "Filename of the parser file",
            "type": str,
            "function_argumemnt": {"pfilename": "pfilename"},
        },
    ]

    generate_sequences_arguments: dict = {
        "ao": "add_original",
        "amto": "add_mutation_to_original",
        "pfilename": "filename",
        "wc": "write_chromosme",
    }

    vcf_reader: FastaReader = None

    def __init__(self, vcf_path: str, fasta_path: str):
        self.vcf_reader = FastaReader(vcf_path, fasta_path)
        logging.info("Loading finalized\n")

    @property
    @abstractmethod
    def name(self) -> str:
        pass

    def get_vcf(self) -> VcfReader:
        """Returns vcf file

        Returns
        -------
        str
            vcf file
        """
        return self.vcf_reader.get_vcf()

    def get_vcf_reader(self) -> FastaReader:
        """Returns vcf reader

        Returns
        -------
        FastaReader
            vcf reader
        """
        return self.vcf_reader

    @staticmethod
    @abstractmethod
    def retrive_sequence(sequence: str) -> Union[tuple, list]:
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
    def retrive_string_sequence(sequence: str) -> str:
        """Gets a string sequence and returns the sequence in a string type

        Parameters
        ----------
        sequence: str
            Sequence in string format

        Returns
        -------
        The sequence in a string format
        """
        pass

    @abstractmethod
    def method(self, sequence: Union[tuple, list], mutation: str) -> Union[tuple, list]:
        """Parse the sequence with a given mutation to a new sequence with different
        noation, for example:
            from ("a", "b", "c") with mutation = "dd" to ("s", "mm", "s")

        Parameters
        ----------
        sequence: tuple, list
            Sequence in a tuple shape
        muatation: str
            Mutation of the nucleotide

        Returns
        -------
        The parsed sequence
        """
        pass

    def original_sequence_to_string(
        self, prefix: str, sequence: list, mutation: str = None
    ) -> str:
        """Generates a string of a sequence by append the prefix, infix and suffix and
        append at the starts a given prefix. If is needed to as the mutation, it can be
        given using the parameter mutation.

        Parameters
        ----------
        prefix : str
            Prefix to append before the original sequence
        sequence: list
            Sequence
        mutation: str
            Mutation of the sequence

        Returns
        -------
        The sequence in a string shape with a given prefix
        """
        if mutation:
            sequence[1] = mutation
        return f'{prefix}{"".join(sequence)}\n'

    @property
    def default_filename(self) -> str:
        """Returns a default filename for the results

        Returns
        -------
        The results default filename
        """
        return f"parsed_{self.name}_data.pvcf"

    @abstractmethod
    def sequence_to_string(
        self,
        original_sequence: str,
        prefix: str,
        sequence: Union[tuple, list],
        mutation: str,
    ) -> str:
        """Creates a string from a sequence and a mutation from the original sequence
        (in a string shape) and a given prefix.

        Parameters
        ----------
        original_sequence : str
            Original sequence in a string shape
        prefix : str
            Prefix to append before the parsed sequence
        sequence: list, tuple
            Sequence
        mutation: str
            Mutation

        Returns
        -------
        The sequence and original sequence with the prefix
        """
        pass

    def generate_sequences(
        self,
        path: str,
        filename: str = False,
        write_chromosme: bool = False,
        add_original: bool = True,
        prefix_length: int = 5,
        suffix_length: int = 5,
        add_mutation_to_original: bool = True,
    ):
        """Generates a file with the sequences mutated using a defined method for parse
        the sequence and the mutation.

        Parameters
        ----------
        path : str
            Path to store the data
        filename : str = default_filename
            Filename of the result file
        write_chromosme : bool = False
            If true adds the chromosme where the sequences being into the file
        add_original : bool = False
            If true adds the original sequence into the file
        """
        if not filename:
            filename = self.default_filename

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
                    mutation = None
                    if add_mutation_to_original:
                        mutation = i.ALT[0].sequence
                    original_sequence = self.original_sequence_to_string(
                        prefix, sequence, mutation=mutation
                    )

                parsed_data_file.write(
                    self.sequence_to_string(
                        original_sequence, prefix, sequence, i.ALT[0].sequence
                    )
                )
