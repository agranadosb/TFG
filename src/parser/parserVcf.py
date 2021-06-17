# -*- coding: utf-8 -*-

import logging
from abc import ABC, abstractmethod
from typing import Union

from src.argumentParser.abstractArguments import AbstractParserArguments
from src.fasta.fastaReader import FastaReader
from src.logging.tqdmLoggingHandler import TqdmLoggingHandler
from tqdm import tqdm
from vcf import Reader as VcfReader


class ParserVcf(AbstractParserArguments, ABC):
    """Parses data from the VCF and FASTA files and prepares that data for a machine
    learning model.

    The parser needs to define a method named `method`, which generates a new sequence
    from the original one.

    For example, if we have the sequence `abbc` (a, b, c) and the mutation `aa` this
    `method` changes it into other representation, for example to `smms` (s, mm, s)
    where the mutation will be added (and transformed) in the middle of a tuple, being
    the start and the end of that tuple the prefix and the suffix respectively:

    ```python
        ("s", "mm", "s") = method(("a", "b", "c"), "bb")
    ```

    The parser will have static methods to send and retrieve the parsed and transformed
    sequence to a file in a specified format. To do this, the parser will need to be defined the methods:

    - ** sequence_to_string **: Transforms a sequence (in a tuple format) to a string
    format, adding to this string the original sequence (in a string format), the
    mutation and a prefix.
    - ** retrive_sequence **: Gets a sequence in a string format and returns it in a
    tuple format. This method will get the string sequences that had been transformed by
    `sequence_to_string` method.
    - ** retrive_string_sequence **: Gets a sequence in a string format and returns it
    in a different string format. This method will get the string sequences that had
    been transformed by `sequence_to_string` method.

    Parameters
    ----------
    vcf_path: str
        Path of the vcf file.
    fasta_path: str
        Path of the fasta file.
    """

    _arguments: list = [
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

    _generate_sequences_arguments: dict = {
        "ao": "add_original",
        "amto": "add_mutation_to_original",
        "pfilename": "filename",
        "wc": "write_chromosme",
        "parser_prefix": "prefix_length",
        "parser_suffix": "suffix_length",
        "result_folder": "path",
    }
    """ Mapping between command line arguments and function arguments of the
    **generate_sequence** method """

    _fasta_reader: FastaReader = None
    """ Fasta reader """

    def __init__(self, vcf_path: str, fasta_path: str):
        logging.info("Loading vcf file")
        self._vcf_file = VcfReader(open(vcf_path, "r"))

        self._fasta_reader = FastaReader(fasta_path)
        logging.info("Loading finalized\n")

    @property
    @abstractmethod
    def name(self) -> str:
        """Name of the parser."""
        pass

    @property
    def _default_filename(self) -> str:
        """Default filename for the results."""
        return f"parsed_{self.name}_data.pvcf"

    def get_vcf(self) -> VcfReader:
        """Returns vcf file.

        Returns
        -------
        Vcf file.
        """
        return self._vcf_file

    def get_fasta_reader(self) -> FastaReader:
        """Returns vcf reader.

        Returns
        -------
        Vcf reader.
        """
        return self._fasta_reader

    def _original_sequence_to_string(
        self,
        prefix: str,
        sequence: list,
        mutation: str = None,
        original_separator: str = "|",
    ) -> str:
        """Generates a string of a sequence by appending the prefix, infix, and suffix
        from a sequence (in a tuple type) and append at the start of the result a given
        prefix.

        If a mutation is given, the mutation is added in the part of the infix and the
        infix is added as a suffix that starts with a given character ("|" by default).

        For example, if we have the sequence **(ACGT,ACGT,ACGT)**, the prefix
        **prefix_name** and the mutation **TGCA**, the method will return:

        ```python
            "prefix_nameACGTTGCAACGT|ACGT"
        ```

        Parameters
        ----------
        prefix: str
            Prefix to append before the original sequence.
        sequence: list
            Sequence.
        mutation: str
            Mutation of the sequence.
        original_separator: str
            Symbol to append to the infix if a mutation is given.

        Returns
        -------
        The sequence in a string shape.
        """
        suffix = ""
        if mutation:
            suffix = f"{original_separator}{sequence[1]}"
            sequence[1] = mutation
        return f'{prefix}{"".join(sequence)}{suffix}\n'

    @abstractmethod
    def method(self, sequence: Union[tuple, list], mutation: str) -> Union[tuple, list]:
        """Parse the sequence with a given mutation to a new sequence with different
        noation.

        For example, if we have the sequence `abbc` (a, b, c) and the mutation `aa` this
        `method` changes it into other representation, for example to `smms` (s, mm, s)
        where the mutation will be added (and transformed) in the middle of a tuple, being
        the start and the end of that tuple the prefix and the suffix respectively:

        ```python
            ("s", "mm", "s") = method(("a", "b", "c"), "bb")
        ```

        Parameters
        ----------
        sequence: tuple, list
            Sequence in a tuple shape.
        muatation: str
            Mutation.

        Returns
        -------
        Transformed sequence.
        """
        pass

    @abstractmethod
    def sequence_to_string(
        self,
        sequence: Union[tuple, list],
        mutation: str,
        original_sequence: str,
        prefix: str,
    ) -> str:
        """Transforms a sequence (in a tuple format) to a string format, adding to this
        string the original sequence (in a string format), the mutation and a prefix.

        Parameters
        ----------
        sequence: list, tuple
            Sequence
        mutation: str
            Mutation
        original_sequence: str
            Original sequence in a string shape
        prefix: str
            Prefix to append before the parsed sequence

        Returns
        -------
        The sequence and original sequence with the prefix
        """
        pass

    def _generate_sequences(
        self,
        path: str,
        filename: str = False,
        write_chromosme: bool = False,
        add_original: bool = True,
        prefix_length: int = 5,
        suffix_length: int = 5,
        add_mutation_to_original: bool = True,
    ):
        """Generates a file with the mutated sequences using the method `method` for
        parse the sequence and the mutation.

        There is many options for add more information bout the sequences, as add the
        source chromsome of the sequence, the length of the prefix and suffix, add the
        original sequence or add the mutation in the original sequence.

        Parameters
        ----------
        path: str
            Path to store the data.
        filename: str = default_filename
            Filename of the result file.
        write_chromosme: bool = False
            If true add the chromosome where the sequences is from into the file.
        add_original: bool = False
            If true adds the original sequence into the file.
        prefix_length: int = 5
            Length of the prefix.
        suffix_length: int = 5
            Length of the suffix.
        add_mutation_to_original: bool = True
            If true Add mutation to original sequence.

        Returns
        -------
        Parsed sequences from vcf.
        """
        logger = logging.getLogger()
        tqdm_out = TqdmLoggingHandler(logger, level=logging.INFO)

        if not filename:
            filename = self._default_filename

        logging.info(f"Parsing sequences using {self.name}")
        sequences = []
        with open(f"{path}/{filename}", "w") as parsed_data_file:
            for i in tqdm(self.get_vcf(), file=tqdm_out):
                sequence = self.get_fasta_reader().get_sequence(
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
                    original_sequence = self._original_sequence_to_string(
                        prefix, sequence.copy(), mutation=mutation
                    )

                parsed_sequence = self.sequence_to_string(
                    sequence, i.ALT[0].sequence, original_sequence, prefix
                )

                sequences.append(parsed_sequence)
                parsed_data_file.write(parsed_sequence)

        logging.info("Parsing finalized\n")
        return sequences

    def generate_sequences(self, **kwargs):
        """Generates a file with the mutated sequences using the method `method` for
        parse the sequence and the mutation.

        There is many options for add more information bout the sequences, as add the
        source chromsome of the sequence, the length of the prefix and suffix, add the
        original sequence or add the mutation in the original sequence.

        Parameters
        ----------
        path: str
            Path to store the data.
        filename: str = default_filename
            Filename of the result file.
        write_chromosme: bool = False
            If true add the chromosome where the sequences is from into the file.
        add_original: bool = False
            If true adds the original sequence into the file.
        prefix_length: int = 5
            Length of the prefix.
        suffix_length: int = 5
            Length of the suffix.
        add_mutation_to_original: bool = True
            If true Add mutation to original sequence.

        Returns
        -------
        Parsed sequences from vcf.
        """
        self._generate_sequences(
            **self.get_generate_sequences_arguments(**kwargs),
        )

    @staticmethod
    @abstractmethod
    def retrive_sequence(sequence: str) -> Union[tuple, list]:
        """Gets a sequence in a string format and returns it in a tuple format. This
        method will get the string sequences that had been transformed by
        `sequence_to_string` method.

        Parameters
        ----------
        sequence: str
            Sequence in string format.

        Returns
        -------
        Sequence in a tuple format.
        """
        pass

    @staticmethod
    @abstractmethod
    def retrive_string_sequence(sequence: str) -> str:
        """Gets a sequence in a string format and returns it in a different string
        format. This method will get the string sequences that had been transformed by
        `sequence_to_string` method.

        Parameters
        ----------
        sequence: str
            Sequence in string format.

        Returns
        -------
        The sequence in a different string format.
        """
        pass
