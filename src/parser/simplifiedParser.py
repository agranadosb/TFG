# -*- coding: utf-8 -*-

from src.parser.parserVcf import ParserVcf


class SimplifiedParserVcf(ParserVcf):
    mutations_symbols = ['m']

    name = "simplified"

    def method(self, sequence: tuple, mutation):
        """Generates the simplified sequence by a equence:
            sequence:               (ACGTGGT,CAA,GTCC)
            sequence_simplified:    (lllllll,mmm,rrrr)

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

        return (
            "l" * len(sequence[0]),
            "m" * len(mutation[0].sequence),
            "r" * len(sequence[2]),
        )
