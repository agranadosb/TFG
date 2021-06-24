# -*- coding: utf-8 -*-

import pathlib
from unittest import TestCase

from src.fasta.fastaReader import FastaReader


class TestFastaReader(TestCase):
    def setUp(self) -> None:
        self.static_dir = f"{str(pathlib.Path(__file__).parent.absolute())}/static/"
        self.reader = FastaReader(
            f"{self.static_dir}test.fa.gz",
        )

        return super().setUp()

    def test__set_fasta_file(self):
        filename = f"{self.static_dir}test.fa"

        result = self.reader.fasta_file.name

        self.assertEqual(result, filename)

    def test_get_sequence(self):
        chromosome = "chr3"
        nucleotide = "GC"
        pos = 9
        from_nuc = 9
        to_nuc = 33
        sequence = ["AGCTAGCTA", "GC", "TAGCTAGCTAGCTAGCTAGCTA"]

        result = self.reader.sequence(chromosome, nucleotide, pos, from_nuc, to_nuc)

        self.assertEqual(result, sequence)


class TestFastaReader1Line(TestCase):
    def setUp(self) -> None:
        self.static_dir = f"{str(pathlib.Path(__file__).parent.absolute())}/static/"
        self.reader = FastaReader(
            f"{self.static_dir}test2.fa.gz",
        )

        return super().setUp()

    def test__set_fasta_file(self):
        filename = f"{self.static_dir}test2.fa"

        result = self.reader.fasta_file.name

        self.assertEqual(result, filename)

    def test__parse_chromosomes(self):
        data = {
            "1": {
                "name": "1",
                "line_length": 29,
                "label_length": 8,
                "index_start": 0,
                "length": 29,
                "labels": ["chr1", "1"],
            },
            "2": {
                "name": "2",
                "line_length": 36,
                "label_length": 8,
                "index_start": 38,
                "length": 36,
                "labels": ["chr2", "2"],
            },
            "3": {
                "name": "3",
                "line_length": 33,
                "label_length": 8,
                "index_start": 83,
                "length": 33,
                "labels": ["chr3", "3"],
            },
        }

        result = self.reader._parse_chromosomes()

        self.assertEqual(result, data)

    def test_get_sequence(self):
        chromosome = "chr3"
        nucleotide = "GC"
        pos = 9
        from_nuc = 9
        to_nuc = 33
        sequence = ["AGCTAGCTA", "GC", "TAGCTAGCTAGCTAGCTAGCTA"]

        result = self.reader.sequence(chromosome, nucleotide, pos, from_nuc, to_nuc)

        self.assertEqual(result, sequence)
