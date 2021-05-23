# -*- coding: utf-8 -*-

from typing import Union

from src.parser.parserVcf import ParserVcf
from src.utils.genomics import generate_dict_values

PREFIX_SYMBOLS = ["q", "w", "e", "r"]
MUTATIONS_SYMBOLS = ["a", "s", "d", "f"]
SUFFIX_SYMBOLS = ["z", "x", "c", "v"]


class ExtendedParserVcf(ParserVcf):
    """Parses data from vcf and fasta and prepare that data for a machine learning model

    The extended parser gets a sequence (prefix, infix, suffix) and maps each character
    to another depending if the character is in the prefix, infix or suffix:

     - Prefix mapping:
         - A -> q
         - C -> w
         - G -> e
         - T -> r
     - Infix mapping:
         - A -> a
         - C -> s
         - G -> d
         - T -> f
    - Suffix mapping:
         - A -> z
         - C -> x
         - G -> c
         - T -> v

    So, for a given sequence ("ACGT", "ACGT", "ACGT") the parsers changes it to:
     - ("qwer", "asdf", "zxcv")

    Parameters
    ----------
    vcf_path : str
        Path of the vcf file
    fasta_path : str
        Path of the fasta file
    """

    prefix_map: dict = generate_dict_values(PREFIX_SYMBOLS)
    mutations_map: dict = generate_dict_values(MUTATIONS_SYMBOLS)
    suffix_map: dict = generate_dict_values(SUFFIX_SYMBOLS)

    mutations_symbols: Union[list, tuple] = MUTATIONS_SYMBOLS

    name: str = "extended"

    def sequence_to_string(
        self,
        original_sequence: str,
        prefix: str,
        sequence: Union[tuple, list],
        mutation: str,
    ) -> str:
        """Creates a string from a sequence and a mutation with the original sequence
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
        Sequence and original sequence with the prefix
        """
        result_sequence = self.method(sequence, mutation)

        parsed_sequence_prefix = "-".join(result_sequence[0])
        parsed_sequence_infix = "-".join(result_sequence[1])
        parsed_sequence_suffix = "-".join(result_sequence[2])

        parsed_sequence = (
            parsed_sequence_prefix,
            parsed_sequence_infix,
            parsed_sequence_suffix,
        )

        return f'{original_sequence}{prefix}{" ".join(parsed_sequence)}\n'

    def method(self, sequence: Union[tuple, list], mutation: str) -> tuple:
        """Generates the extended sequence from a equence:
            sequence:               (ACGTGGT,CAA,GTCC)
            sequence extended:      ([a,s,d,f,d,d,f],[e,q,q],[c,v,x,x])

        Parameters
        ----------
        sequence : tuple
            Sequence to simplify
        mutation : str
            Mutation sequence

        Returns
        -------
        tuple
            sequence simplified
        """

        left = [self.prefix_map[nucletid.upper()] for nucletid in sequence[0]]
        middle = [self.mutations_map[nucletid.upper()] for nucletid in mutation]
        right = [self.suffix_map[nucletid.upper()] for nucletid in sequence[2]]

        return (left, middle, right)

    @staticmethod
    def retrive_sequence(sequence: str) -> tuple:
        """Gets a string sequence and returns the sequence in a tuple type

        Parameters
        ----------
        sequence: str
            Sequence in string format

        Returns
        -------
        Sequence in a list format
        """
        sequence_diveded = sequence.split(" ")

        prefix = sequence_diveded[0].split("-")
        infix = sequence_diveded[1].split("-")
        suffix = sequence_diveded[2].split("-")

        return (
            "".join(prefix).rstrip(),
            "".join(infix).rstrip(),
            "".join(suffix).rstrip(),
        )

    @staticmethod
    def retrive_string_sequence(sequence: str) -> str:
        """Gets a string sequence and returns the sequence in a tuple type

        Parameters
        ----------
        sequence: str
            Sequence in string format

        Returns
        -------
        Sequence in a string format
        """
        sequence_diveded = sequence.split(" ")

        prefix = sequence_diveded[0].split("-")
        infix = sequence_diveded[1].split("-")
        suffix = sequence_diveded[2].split("-")

        return f"""{"".join(prefix).rstrip()}{"".join(infix).rstrip()}{"".join(suffix).rstrip()}"""
