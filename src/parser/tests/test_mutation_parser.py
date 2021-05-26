# -*- coding: utf-8 -*-

import pathlib
from unittest import TestCase

from src.parser.mutationParser import MutationParser


class TestMutationParser(TestCase):
    def setUp(self) -> None:
        self.static_dir = f"{str(pathlib.Path(__file__).parent.absolute())}/static/"
        self.parser = MutationParser(
            f"{self.static_dir}vcfTest.vcf",
            f"{self.static_dir}test.fa.gz",
        )

        return super().setUp()

    def test_method_insertion(self):
        sequence = ("ACGT", "", "ACGT")
        mutation = "ACGT"
        mutation_type_sequence = (
            ["q", "w", "e", "r"],
            ["t", "y", "u", "i"],
            ["z", "x", "c", "v"],
        )

        result = self.parser.method(sequence, mutation)

        self.assertEqual(result, mutation_type_sequence)

    def test_method_delete(self):
        sequence = ("ACGT", "ACGT", "ACGT")
        mutation = ""
        mutation_type_sequence = (
            ["q", "w", "e", "r"],
            ["g", "h", "j", "k"],
            ["z", "x", "c", "v"],
        )

        result = self.parser.method(sequence, mutation)

        self.assertEqual(result, mutation_type_sequence)

    def test_method_replace(self):
        sequence = ("ACGT", "ACGT", "ACGT")
        mutation = "TGCA"
        mutation_type_sequence = (
            ["q", "w", "e", "r"],
            ["a", "s", "d", "f"],
            ["z", "x", "c", "v"],
        )

        result = self.parser.method(sequence, mutation)

        self.assertEqual(result, mutation_type_sequence)

    def test_method_mix(self):
        sequence = ("ACGT", "AGCT", "ACGT")
        mutation = "GTTCAC"
        mutation_type_sequence = (
            ["q", "w", "e", "r"],
            ["u", "a", "d", "t", "f"],
            ["z", "x", "c", "v"],
        )

        result = self.parser.method(sequence, mutation)

        self.assertEqual(result, mutation_type_sequence)
