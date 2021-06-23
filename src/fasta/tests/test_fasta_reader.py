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

    def test__chromosome_index(self):
        chromosomes = [("chr1", 0, 12), ("chr2", 4, 12), ("chr3", 8, 12)]
        indexes = [0, 38, 83]

        result = [
            self.reader._chromosome_index(*chromosome) for chromosome in chromosomes
        ]

        self.assertEqual(result, indexes)

    def test__parse_chromosomes(self):
        data = {
            "chr1": {
                "label_length": 6,
                "index_start": 0,
                "length": 29,
            },
            "chr2": {
                "label_length": 6,
                "index_start": 38,
                "length": 36,
            },
            "chr3": {
                "label_length": 6,
                "index_start": 83,
                "length": 33,
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
