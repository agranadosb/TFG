# -*- coding: utf-8 -*-

import pathlib
from unittest import TestCase

from src.parser.extendedParser import ExtendedParserVcf
from src.parser.tests.factories import SequenceFaker


class TestExtendedParserVcf(TestCase):
    def setUp(self) -> None:
        self.static_dir = f"{str(pathlib.Path(__file__).parent.absolute())}/static/"
        self.parser = ExtendedParserVcf(
            f"{self.static_dir}vcfTest.vcf",
            f"{self.static_dir}test.fa.gz",
        )

        return super().setUp()

    def test_method(self):
        sequence = ("ACGTGGT", "CAA", "GTCC")
        mutation = SequenceFaker("AAA")
        prefix = ["01", "02", "03", "04", "03", "03", "04"]
        infix = ["11", "11", "11"]
        suffix = ["23", "24", "22", "22"]

        result = self.parser.method(sequence, [mutation])

        self.assertEqual(result[0], prefix)
        self.assertEqual(result[1], infix)
        self.assertEqual(result[2], suffix)
