# -*- coding: utf-8 -*-

import pathlib
from unittest import TestCase

from src.parser.extendedParser import ExtendedParserVcf
from src.parser.tests.factories import SequenceFaker
from src.parser.tests.helpers import random_string


class TestExtendedParserVcf(TestCase):
    def setUp(self) -> None:
        self.static_dir = f"{str(pathlib.Path(__file__).parent.absolute())}/static/"
        self.parser = ExtendedParserVcf(
            f"{self.static_dir}vcfTest.vcf",
            f"{self.static_dir}test.fa.gz",
        )

        return super().setUp()

    def test_sequence_to_string(self):
        original_sequence = random_string()
        prefix = random_string()
        sequence = ("ACGTGGT", "CAA", "GTCC")
        sequence_prefix = "a-s-d-f-d-d-f"
        sequence_infix = "q-q-q"
        sequence_suffix = "c-v-x-x"
        sequence_result = f"{sequence_prefix} {sequence_infix} {sequence_suffix}"
        mutation = SequenceFaker("AAA")
        sequence_string = f"{original_sequence}{prefix}{sequence_result}\n"

        result = self.parser.sequence_to_string(
            original_sequence, prefix, sequence, [mutation]
        )

        self.assertEqual(result, sequence_string)

    def test_retrive_string_sequence(self):
        sequence = ("ACGTGGT", "CAA", "GTCC")
        sequence_prefix = "a-s-d-f-d-d-f"
        sequence_infix = "q-q-q"
        sequence_suffix = "c-v-x-x"
        sequence = f"{sequence_prefix} {sequence_infix} {sequence_suffix}"
        result_prefix = "asdfddf"
        result_infix = "qqq"
        result_suffix = "cvxx"
        result_sequence = f"{result_prefix}{result_infix}{result_suffix}"

        result = ExtendedParserVcf.retrive_string_sequence(sequence)

        self.assertEqual(result, result_sequence)

    def test_method(self):
        sequence = ("ACGTGGT", "CAA", "GTCC")
        mutation = SequenceFaker("AAA")
        prefix = ["a", "s", "d", "f", "d", "d", "f"]
        infix = ["q", "q", "q"]
        suffix = ["c", "v", "x", "x"]

        result = self.parser.method(sequence, [mutation])

        self.assertEqual(result[0], prefix)
        self.assertEqual(result[1], infix)
        self.assertEqual(result[2], suffix)
