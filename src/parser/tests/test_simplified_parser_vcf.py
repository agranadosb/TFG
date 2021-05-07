# -*- coding: utf-8 -*-

import pathlib
from unittest import TestCase

from src.parser.simplifiedParser import SimplifiedParserVcf
from src.parser.tests.factories import SequenceFaker
from src.parser.tests.helpers import random_string


class TestSimplifiedParserVcf(TestCase):
    def setUp(self) -> None:
        self.static_dir = f"{str(pathlib.Path(__file__).parent.absolute())}/static/"
        self.parser = SimplifiedParserVcf(
            f"{self.static_dir}vcfTest.vcf",
            f"{self.static_dir}test.fa.gz",
        )

        return super().setUp()

    def test_sequence_to_string(self):
        original_sequence = random_string()
        prefix = random_string()
        sequence = ("ACGTGGT", "CAA", "GTCC")
        sequence_result = "lllllllmmmrrrr"
        mutation = SequenceFaker("AAA")
        sequence_string = f"{original_sequence}{prefix}{sequence_result}\n"

        result = self.parser.sequence_to_string(
            original_sequence, prefix, sequence, [mutation]
        )

        self.assertEqual(result, sequence_string)

    def test_method(self):
        sequence = ("ACGTGGT", "CAA", "GTCC")
        mutation = SequenceFaker("AAA")
        prefix = "lllllll"
        infix = "mmm"
        suffix = "rrrr"

        result = self.parser.method(sequence, [mutation])

        self.assertEqual(result[0], prefix)
        self.assertEqual(result[1], infix)
        self.assertEqual(result[2], suffix)
