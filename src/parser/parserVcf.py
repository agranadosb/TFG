# -*- coding: utf-8 -*-

import logging
from abc import ABC, abstractmethod
from typing import Union

from src.argumentParser.abstractArguments import AbstractParserArguments
from src.logging.tqdmLoggingHandler import TqdmLoggingHandler
from src.fasta.fastaReader import FastaReader
from tqdm import tqdm
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
     - **retrive_sequence**: Get a sequence in a string shape and returns it as a tuple
     - **retrive_string_sequence**: Get a sequence in a string shape and returns it as a
    string

    Parameters
    ----------
    vcf_path: str
        Path of the vcf file
    fasta_path: str
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
    """ Arguments that will be used by command line """

    generate_sequences_arguments: dict = {
        "ao": "add_original",
        "amto": "add_mutation_to_original",
        "pfilename": "filename",
        "wc": "write_chromosme",
    }
    """ Mapping between command line arguments and function arguments of the
    **generate_sequence** method """

    fasta_reader: FastaReader = None
    """ Fasta reader """

    def __init__(self, vcf_path: str, fasta_path: str):
        logging.info("Loading vcf file")
        self.vcf_file = VcfReader(open(vcf_path, "r"))

        self.fasta_reader = FastaReader(fasta_path)
        logging.info("Loading finalized\n")

    @property
    @abstractmethod
    def name(self) -> str:
        """Name of the parser"""
        pass

    @property
    def default_filename(self) -> str:
        """Default filename for the results"""
        return f"parsed_{self.name}_data.pvcf"

    def get_vcf(self) -> VcfReader:
        """Returns vcf file

        Returns
        -------
        Vcf file
        """
        return self.vcf_file

    def get_fasta_reader(self) -> FastaReader:
        """Returns vcf reader

        Returns
        -------
        Vcf reader
        """
        return self.fasta_reader

    def original_sequence_to_string(
        self,
        prefix: str,
        sequence: list,
        mutation: str = None,
        original_separator: str = "|",
    ) -> str:
        """Generates a string of a sequence by append the prefix, infix and suffix and
        append at the starts a given prefix. If is needed to as the mutation, it can be
        given using the parameter mutation.

        For example, if we have the sequence **(ACGT,ACGT,ACGT)**, the prefix **prefix**
        and the mutation **TGCA**, the method will return:

            "prefixACGTTGCAACGT|ACGT"

        Parameters
        ----------
        prefix: str
            Prefix to append before the original sequence
        sequence: list
            Sequence
        mutation: str
            Mutation of the sequence
        original_separator: str
            Symbol that divides an original sequence between the sequence and the
            reference sequence

        Returns
        -------
        The sequence in a string shape with a given prefix
        """
        suffix = ""
        if mutation:
            suffix = f"{original_separator}{sequence[1]}"
            sequence[1] = mutation
        return f'{prefix}{"".join(sequence)}{suffix}\n'

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
        original_sequence: str
            Original sequence in a string shape
        prefix: str
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
        path: str
            Path to store the data
        filename: str = default_filename
            Filename of the result file
        write_chromosme: bool = False
            If true adds the chromosme where the sequences being into the file
        add_original: bool = False
            If true adds the original sequence into the file
        prefix_length: int = 5
            Length of the prefix
        suffix_length: int = 5
            Length of the suffix
        add_mutation_to_original: bool = True
            Add mutation to original sequence if true
        
        Returns
        -------
        Parsed sequences from vcf
        """
        logger = logging.getLogger()
        tqdm_out = TqdmLoggingHandler(logger, level=logging.INFO)

        if not filename:
            filename = self.default_filename

        logging.info(f"Parsing sequences using {self.name}")
        sequences = []
        with open(f"{path}/{filename}", "w") as parsed_data_file:
            for i in tqdm(self.get_vcf(), file=tqdm_out):
                sequence = self.fasta_reader.get_sequence(
                    i.CHROM, i.REF, i.POS - 1, prefix_length, suffix_length
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
                        prefix, sequence.copy(), mutation=mutation
                    )

                parsed_sequence = self.sequence_to_string(
                    original_sequence, prefix, sequence, i.ALT[0].sequence
                )

                sequences.append(parsed_sequence)
                parsed_data_file.write(parsed_sequence)

        logging.info("Parsing finalized\n")
        return sequences

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
