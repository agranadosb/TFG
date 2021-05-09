# -*- coding: utf-8 -*-

from src.parser.parserVcf import ParserVcf

NUCLETODIES = ["A", "C", "G", "T"]


def generate_dict_values(function):
    return {key: function(value) for (key, value) in zip(NUCLETODIES, range(4))}


class ExtendedParserVcf(ParserVcf):

    left_map = generate_dict_values(lambda value: f"0{str(value + 1)}")
    middle_map = generate_dict_values(lambda value: str(value + 11))
    right_map = generate_dict_values(lambda value: str(value + 21))

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
            sequence_simplified:    ([01,02,03,04,03,03,04],[12,11,11],[23,24,22,22])

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
