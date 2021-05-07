# -*- coding: utf-8 -*-

from src.parser.parserVcf import ParserVcf

NUCLETODIES = ["A", "C", "G", "T"]


class ExtendedParserVcf(ParserVcf):

    left_map = {key: value + 1 for (key, value) in zip(NUCLETODIES, range(4))}
    middle_map = {key: value + 10 + 1 for (key, value) in zip(NUCLETODIES, range(4))}
    right_map = {key: value + 20 + 1 for (key, value) in zip(NUCLETODIES, range(4))}

    name = "extended"

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

        left = [self.left_map[nucletid] for nucletid in sequence[0]]
        middle = [self.middle_map[nucletid] for nucletid in mutation[0].sequence]
        right = [self.right_map[nucletid] for nucletid in sequence[2]]

        return (left, middle, right)
