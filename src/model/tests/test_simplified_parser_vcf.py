# -*- coding: utf-8 -*-

import pathlib
from unittest import TestCase

from src.model.simplifiedParser import SimplifiedParserVcf
from src.model.tests.factories import SequenceFaker


class TestSimplifiedParserVcf(TestCase):
    def setUp(self) -> None:
        self.static_dir = f"{str(pathlib.Path(__file__).parent.absolute())}/static/"
        self.parser = SimplifiedParserVcf(
            f"{self.static_dir}vcfTest.vcf",
            f"{self.static_dir}test.fa.gz",
        )

        return super().setUp()

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
