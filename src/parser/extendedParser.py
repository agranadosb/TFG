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
    vcf_path: str
        Path of the vcf file
    fasta_path: str
        Path of the fasta file
    """

    prefix_map: dict = generate_dict_values(PREFIX_SYMBOLS)
    """ Mapping between nucleotides and the prefix symbols """

    mutations_map: dict = generate_dict_values(MUTATIONS_SYMBOLS)
    """ Mapping between nucleotides and the mutation symbols """

    suffix_map: dict = generate_dict_values(SUFFIX_SYMBOLS)
    """ Mapping between nucleotides and the suffix symbols """

    mutations_symbols: Union[list, tuple] = MUTATIONS_SYMBOLS
    """ List of mutation symbols """

    @classmethod
    def inverse_mutations_map(cls):
        """ Inverse of mutations_map """
        return {value: key for key, value in cls.mutations_map.items()}

    name: str = "extended"

    @classmethod
    def method(cls, sequence: Union[tuple, list], mutation: str) -> tuple:
        """Generates the extended sequence from a equence:
            
            sequence:               (ACGT,ACGT,ACGT)
            sequence extended:      ([q,w,e,r],[a,s,d,f],[z,x,c,v])

        Parameters
        ----------
        sequence: tuple
            Sequence to simplify
        mutation: str
            Mutation sequence

        Returns
        -------
        Sequence simplified
        """

        left = [cls.prefix_map[nucletid.upper()] for nucletid in sequence[0]]
        middle = [cls.mutations_map[nucletid.upper()] for nucletid in mutation]
        right = [cls.suffix_map[nucletid.upper()] for nucletid in sequence[2]]

        return (left, middle, right)

    def sequence_to_string(
        self,
        original_sequence: str,
        prefix: str,
        sequence: Union[tuple, list],
        mutation: str,
        separator_symbols: str = "-",
        separator_sequences: str = " ",
        prefix_separator: str = "*",
    ) -> str:
        """Creates a string from a sequence and a mutation with the original sequence
        (in a string shape) and a given prefix.

        For example, if our sequence is **(ACGT, ACGT, ACGT)**, te original sequence is
        **ACGT**, the prefix is **prefix** and the mutation is **TGCA**, the method will
        return:

            "ACGT*prefix*q-w-e-r f-d-s-a z-x-c-v"

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
        separator_symbols: str = "-"
            Symbol that is between each symbol of the parsed sequence
        separator_sequences: str = " "
            Symbol that is between each part of the sequence (prefix, symbol, infix,
            symbol, suffix)
        prefix_separator: str = "*"
            Symbols that divide each part of the prefix

        Returns
        -------
        Sequence and original sequence with the prefix
        """
        result_sequence = self.method(sequence, mutation)

        parsed_sequence_prefix = separator_symbols.join(result_sequence[0])
        parsed_sequence_infix = separator_symbols.join(result_sequence[1])
        parsed_sequence_suffix = separator_symbols.join(result_sequence[2])

        parsed_sequence = (
            parsed_sequence_prefix,
            parsed_sequence_infix,
            parsed_sequence_suffix,
        )

        string_sequence_prefix = (
            f"{original_sequence}{prefix_separator}{prefix}{prefix_separator}"
        )

        string_sequence = separator_sequences.join(parsed_sequence)

        return f"{string_sequence_prefix}{string_sequence}\n"

    @staticmethod
    def retrive_sequence(
        sequence: str,
        separator_symbols: str = "-",
        separator_sequences: str = " ",
        prefix_separator: str = "*",
    ) -> tuple:
        """Gets a string sequence and returns the sequence in a tuple type.

        For example, if our sequence is **ACGT\*prefix\*q-w-e-r f-d-s-a z-x-c-v** the
        method will return:

            ("qwer", "fdsa", "zxcv")

        Parameters
        ----------
        sequence: str
            Sequence in string format
        separator_symbols: str = "-"
            Symbol that is between each symbol of the parsed sequence
        separator_sequences: str = " "
            Symbol that is between each part of the sequence (prefix, symbol, infix,
            symbol, suffix)
        prefix_separator: str = "*"
            Symbols that divide each part of the prefix

        Returns
        -------
        Sequence in a list format
        """
        prefix_divided = sequence.split(prefix_separator)
        sequence_diveded = prefix_divided[len(prefix_divided) - 1].split(
            separator_sequences
        )

        prefix = sequence_diveded[0].split(separator_symbols)
        infix = sequence_diveded[1].split(separator_symbols)
        suffix = sequence_diveded[2].split(separator_symbols)

        return (
            "".join(prefix).rstrip(),
            "".join(infix).rstrip(),
            "".join(suffix).rstrip(),
        )

    @staticmethod
    def retrive_string_sequence(
        sequence: str,
        separator_symbols: str = "-",
        separator_sequences: str = " ",
        prefix_separator: str = "*",
    ) -> str:
        """Gets a string sequence and returns the sequence in a tuple type.

        For example, if our sequence is **ACGT\*prefix\*q-w-e-r f-d-s-a z-x-c-v** the
        method will return:

            "qwerfdsazxcv"

        Parameters
        ----------
        sequence: str
            Sequence in string format
        separator_symbols: str = "-"
            Symbol that is between each symbol of the parsed sequence
        separator_sequences: str = " "
            Symbol that is between each part of the sequence (prefix, symbol, infix,
            symbol, suffix)
        prefix_separator: str = "*"
            Symbols that divide each part of the prefix

        Returns
        -------
        Sequence in a string format
        """
        prefix_divided = sequence.split(prefix_separator)
        sequence_diveded = prefix_divided[len(prefix_divided) - 1].split(
            separator_sequences
        )

        prefix = sequence_diveded[0].split(separator_symbols)
        infix = sequence_diveded[1].split(separator_symbols)
        suffix = sequence_diveded[2].split(separator_symbols)

        return f"""{"".join(prefix).rstrip()}{"".join(infix).rstrip()}{"".join(suffix).rstrip()}"""
