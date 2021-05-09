# -*- coding: utf-8 -*-

from src.parser.parserVcf import ParserVcf

NUCLETODIES = ["A", "C", "G", "T"]

PREFIX_SYMBOLS = ["a", "s", "d", "f"]
INFIX_SYMBOLS = ["q", "w", "e", "r"]
SUFFIX_SYMBOLS = ["z", "x", "c", "v"]


def generate_dict_values(lst):
    return {key: value for (key, value) in zip(NUCLETODIES, lst)}


class ExtendedParserVcf(ParserVcf):

    left_map = generate_dict_values(PREFIX_SYMBOLS)
    middle_map = generate_dict_values(INFIX_SYMBOLS)
    right_map = generate_dict_values(SUFFIX_SYMBOLS)

    name = "extended"

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

    def method(self, sequence: tuple, mutation: str):
        """Generates the simplified sequence by a equence:
            sequence:               (ACGTGGT,CAA,GTCC)
            sequence_simplified:    ([a,s,d,f,d,d,f],[e,q,q],[c,v,x,x])

        Parameters
        ----------
        sequence : tuple
            Sequence to simplify
        mutation : tuple
            Mutation sequence

        Returns
        -------
        tuple
            sequence simplified
        """

        left = [self.left_map[nucletid.upper()] for nucletid in sequence[0]]
        middle = [
            self.middle_map[nucletid.upper()] for nucletid in mutation[0].sequence
        ]
        right = [self.right_map[nucletid.upper()] for nucletid in sequence[2]]

        return (left, middle, right)

    def retrive_string_sequence(self, sequence):
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

        return ("".join(prefix), "".join(infix), "".join(suffix))
