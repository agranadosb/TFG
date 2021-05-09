# -*- coding: utf-8 -*-

from src.parser.parserVcf import ParserVcf

NUCLETODIES = ["A", "C", "G", "T"]

def generate_dict_values(function):
    return {
        key: function(value) for (key, value) in zip(NUCLETODIES, range(4))
    }

class ExtendedParserVcf(ParserVcf):

    left_map = generate_dict_values(lambda value: f"0{str(value + 1)}")
    middle_map = generate_dict_values(lambda value: str(value + 11))
    right_map = generate_dict_values(lambda value: str(value + 21))

    name = "extended"

    def sequence_to_string(self, original_sequence, prefix, sequence, mutation):

        result_sequence = self.method(sequence, mutation)
        result_sequence_prefix = "".join(result_sequence[0])
        return f'{original_sequence}{prefix}{"".join()}\n'

    def method(self, sequence: tuple, mutation: str):
        """Generates the simplified sequence by a equence:
            sequence:               (ACGTGGT,CAA,GTCC)
            sequence_simplified:    ([1,2,3,4,3,3,4],[12,11,11],[23,24,22,22])

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
        middle = [self.middle_map[nucletid.upper()] for nucletid in mutation[0].sequence]
        right = [self.right_map[nucletid.upper()] for nucletid in sequence[2]]

        return (left, middle, right)
