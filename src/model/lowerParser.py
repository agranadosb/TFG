# -*- coding: utf-8 -*-

from src.model.parserVcf import ParserVcf


class LowerParserVcf(ParserVcf):
    name = 'lower'

    def method(self, sequence: tuple, mutation: str):
        """Generates the same sequence but in a lower case in the mutations part:
            sequence:               (ACGTGGT,CAA,GTCC)
            sequence_lower:         (ACGTGGT,caa,GTCC)

        Parameters
        ----------
        sequence : tuple
            Sequence to simplify
        mutation : tuple
            Mutation sequence

        Returns
        -------
        tuple
            lower sequence
        """
        return (sequence[0], mutation[0].sequence.lower(), sequence[2])
