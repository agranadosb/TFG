# -*- coding: utf-8 -*-

import pathlib
from logging import setLogRecordFactory
from unittest import TestCase

from src.vcf.vcfReader import FastaReader


class TestFastaReader(TestCase):
    def setUp(self) -> None:
        self.static_dir = f"{str(pathlib.Path(__file__).parent.absolute())}/static/"
        self.reader = FastaReader(
            f"{self.static_dir}vcfTest.vcf",
            f"{self.static_dir}test.fa.gz",
        )

        return super().setUp()

    def test_get_vcf(self):
        filename = f"{self.static_dir}vcfTest.vcf"

        result = self.reader.get_vcf().filename

        self.assertEqual(result, filename)

    def test_get_fasta(self):
        filename = f"{self.static_dir}test.fa"

        result = self.reader.get_fasta().name

        self.assertEqual(result, filename)

    def test_set_fasta_line_length(self):
        length = 12

        result = self.reader.set_fasta_line_length()

        self.assertEqual(result, length)

    def test_set_fasta_file(self):
        filename = f"{self.static_dir}test.fa"

        result = self.reader.get_fasta().name

        self.assertEqual(result, filename)

    def test_get_chromosome_index(self):
        chromosmes = ["chr1", "chr2", "chr3"]
        indexes = [0, 38, 83]

        result = [
            self.reader.get_chromosome_index(chromosome) for chromosome in chromosmes
        ]

        self.assertEqual(result, indexes)

    def test_set_chromosme_start_data(self):
        chromosmes = ["1:>chr1", "5:>chr2", "9:>chr3"]
        data = [
            {"label_length": 6, "line_start": 0, "index_start": 0},
            {"label_length": 6, "line_start": 4, "index_start": 38},
            {"label_length": 6, "line_start": 8, "index_start": 83},
        ]

        result = [
            self.reader.set_chromosme_start_data(chromosme) for chromosme in chromosmes
        ]

        self.assertEqual(result, data)

    def test_set_chromosme_end_data(self):
        chromosomes = 3
        data = [
            {
                "label_length": 6,
                "line_start": 0,
                "index_start": 0,
                "index_ends": 38,
                "next_label_length": 6,
                "is_last": False,
                "index": 0,
                "chromosme_length": 29,
            },
            {
                "label_length": 6,
                "line_start": 4,
                "index_start": 38,
                "index_ends": 83,
                "next_label_length": 6,
                "is_last": False,
                "index": 1,
                "chromosme_length": 36,
            },
            {
                "label_length": 6,
                "line_start": 8,
                "index_start": 83,
                "index_ends": 125,
                "next_label_length": 0,
                "is_last": True,
                "index": 2,
                "chromosme_length": 33,
            },
        ]

        result = [
            self.reader.set_chromosme_end_data(index) for index in range(chromosomes)
        ]

        self.assertEqual(result, data)

    def test_get_fasta_data(self):
        data = {
            "chr1": {
                "label_length": 6,
                "line_start": 0,
                "index_start": 0,
                "index_ends": 38,
                "next_label_length": 6,
                "is_last": False,
                "index": 0,
                "chromosme_length": 29,
            },
            "chr2": {
                "label_length": 6,
                "line_start": 4,
                "index_start": 38,
                "index_ends": 83,
                "next_label_length": 6,
                "is_last": False,
                "index": 1,
                "chromosme_length": 36,
            },
            "chr3": {
                "label_length": 6,
                "line_start": 8,
                "index_start": 83,
                "index_ends": 125,
                "next_label_length": 0,
                "is_last": True,
                "index": 2,
                "chromosme_length": 33,
            },
        }

        result = self.reader.get_fasta_data()

        self.assertEqual(result, data)

    def test_get_nucleotide_index_pos_greater_length(self):
        pos = 1000
        chromosome = "chr1"
        error = False

        try:
            self.reader.get_nucleotide_index(pos, chromosome)
        except IndexError:
            error = True

        self.assertTrue(error)

    def test_get_nucleotide_index_pos_lower_0(self):
        pos = -1
        chromosome = "chr1"
        error = False

        try:
            self.reader.get_nucleotide_index(pos, chromosome)
        except IndexError:
            error = True

        self.assertTrue(error)

    def test_get_nucleotide_index(self):
        pos = 24
        chromosome = "chr1"
        index = 32

        result = self.reader.get_nucleotide_index(pos, chromosome)

        self.assertEqual(result, index)

    def test_get_from_interval(self):
        index = 6
        length = 29
        sequence = "GCATGCATGCATGCATGCATGCATGCATG"

        result = self.reader.get_from_interval(index, length)

        self.assertEqual(result, sequence)

    def test_get_prefix_at_start(self):
        pos = 5
        length = 7
        chromosome = "chr1"
        prefix = "GCATG"

        result = self.reader.get_prefix(pos, length, chromosome)

        self.assertEqual(result, prefix)

    def test_get_prefix(self):
        pos = 13
        length = 7
        chromosome = "chr1"
        prefix = "ATGCATG"

        result = self.reader.get_prefix(pos, length, chromosome)

        self.assertEqual(result, prefix)

    def test_get_prefix_at_end(self):
        pos = 33
        length = 7
        chromosome = "chr2"
        suffix = "AC"

        result = self.reader.get_suffix(pos, length, chromosome)

        self.assertEqual(result, suffix)

    def test_get_prefix(self):
        pos = 22
        length = 7
        chromosome = "chr2"
        suffix = "CTGACTG"

        result = self.reader.get_suffix(pos, length, chromosome)

        self.assertEqual(result, suffix)

    def test_get_nucleotides(self):
        chromosome = "chr3"
        pos = 20
        length = 10
        sequence = "AGCTAGCTAG"

        result = self.reader.get_nucleotides(chromosome, pos, length)

        self.assertEqual(result, sequence)

    def test_get_sequence(self):
        chromosome = "chr3"
        nucleotide = "GC"
        pos = 9
        from_nuc = 200
        to_nuc = 200
        sequence = ["AGCTAGCTA", "GC", "TAGCTAGCTAGCTAGCTAGCTA"]

        result = self.reader.get_sequence(chromosome, nucleotide, pos, from_nuc, to_nuc)

        self.assertEqual(result, sequence)
