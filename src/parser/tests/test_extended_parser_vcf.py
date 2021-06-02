# -*- coding: utf-8 -*-

import pathlib
from unittest import TestCase
from unittest.mock import mock_open, patch

from src.parser.extendedParser import ExtendedParserVcf


class TestExtendedParserVcf_ParserVcf(TestCase):
    def setUp(self) -> None:
        self.static_dir = f"{str(pathlib.Path(__file__).parent.absolute())}/static/"
        self.parser = ExtendedParserVcf(
            f"{self.static_dir}vcfTest.vcf",
            f"{self.static_dir}test.fa.gz",
        )

        return super().setUp()

    def test_name(self):
        name = "extended"

        result = self.parser.name

        self.assertEqual(result, name)

    def test_default_filename(self):
        filename = "parsed_extended_data.pvcf"

        result = self.parser.default_filename

        self.assertEqual(result, filename)

    def test_get_vcf(self):
        filename = f"{self.static_dir}vcfTest.vcf"

        result = self.parser.get_vcf().filename

        self.assertEqual(result, filename)

    def test_get_vcf_reader(self):
        filename = f"{self.static_dir}vcfTest.vcf"

        result = self.parser.get_vcf_reader().get_vcf().filename

        self.assertEqual(result, filename)

    def test_original_sequence_to_string(self):
        prefix = "prefix"
        sequence = ["ACGT", "ACGT", "ACGT"]
        string_sequence = f"{prefix}ACGTACGTACGT\n"

        result = self.parser.original_sequence_to_string(prefix, sequence)

        self.assertEqual(result, string_sequence)

    def test_original_sequence_to_string_with_mutation(self):
        prefix = "prefix"
        sequence = ["ACGT", "ACGT", "ACGT"]
        mutation = "TGCA"
        string_sequence = f"{prefix}ACGTTGCAACGT|ACGT\n"

        result = self.parser.original_sequence_to_string(prefix, sequence, mutation)

        self.assertEqual(result, string_sequence)

    def test_method(self):
        sequence = ("ACGT", "ACGT", "ACGT")
        mutation = "TGCA"
        parsed_sequence = (
            ["q", "w", "e", "r"],
            ["f", "d", "s", "a"],
            ["z", "x", "c", "v"],
        )

        result = self.parser.method(sequence, mutation)

        self.assertEqual(result, parsed_sequence)

    def test_sequence_to_string(self):
        original_sequence = "ACGT"
        prefix = "prefix"
        sequence = ("ACGT", "ACGT", "ACGT")
        mutation = "TGCA"
        parsed_sequence = f"{original_sequence}*{prefix}*q-w-e-r f-d-s-a z-x-c-v\n"

        result = self.parser.sequence_to_string(
            original_sequence, prefix, sequence, mutation
        )

        self.assertEqual(result, parsed_sequence)

    def test_generate_sequences(self):
        sequences = [
            "ACATGC|G\n** a x-z-v-c-x\n",
            "GTTGCAT|CA\n**e f v-c-x-z-v\n",
            "ATGCAACATGC|TG\n**q-r-e-w-q a x-z-v-c-x\n",
            "TGCATTACATGC|G\n**r-e-w-q-r f-a x-z-v-c-x\n",
            "ATGCACG|T\n**q-r-e-w-q s c\n",
            "AGACTG|T\n** a c-z-x-v-c\n",
            "TTCTGAC|GA\n**r f x-v-c-z-x\n",
            "ACTGAAGACTG|CT\n**q-w-r-e-q a c-z-x-v-c\n",
            "CTGACTAGACTG|T\n**w-r-e-q-w f-a c-z-x-v-c\n",
            "CTGACCGACTG|T\n**w-r-e-q-w s c-z-x-v-c\n",
        ]
        with patch("builtins.open", mock_open(read_data="data")):
            result = self.parser.generate_sequences("")

        self.assertEqual(result, sequences)

    def test_generate_sequences_write_chromosme_true(self):
        sequences = [
            "chr1\tACATGC|G\n*chr1\t* a x-z-v-c-x\n",
            "chr1\tGTTGCAT|CA\n*chr1\t*e f v-c-x-z-v\n",
            "chr1\tATGCAACATGC|TG\n*chr1\t*q-r-e-w-q a x-z-v-c-x\n",
            "chr1\tTGCATTACATGC|G\n*chr1\t*r-e-w-q-r f-a x-z-v-c-x\n",
            "chr1\tATGCACG|T\n*chr1\t*q-r-e-w-q s c\n",
            "chr2\tAGACTG|T\n*chr2\t* a c-z-x-v-c\n",
            "chr2\tTTCTGAC|GA\n*chr2\t*r f x-v-c-z-x\n",
            "chr2\tACTGAAGACTG|CT\n*chr2\t*q-w-r-e-q a c-z-x-v-c\n",
            "chr2\tCTGACTAGACTG|T\n*chr2\t*w-r-e-q-w f-a c-z-x-v-c\n",
            "chr2\tCTGACCGACTG|T\n*chr2\t*w-r-e-q-w s c-z-x-v-c\n",
        ]
        with patch("builtins.open", mock_open(read_data="data")):
            result = self.parser.generate_sequences("", write_chromosme=True)

        self.assertEqual(result, sequences)

    def test_generate_sequences_add_original_false(self):
        sequences = [
            "** a x-z-v-c-x\n",
            "**e f v-c-x-z-v\n",
            "**q-r-e-w-q a x-z-v-c-x\n",
            "**r-e-w-q-r f-a x-z-v-c-x\n",
            "**q-r-e-w-q s c\n",
            "** a c-z-x-v-c\n",
            "**r f x-v-c-z-x\n",
            "**q-w-r-e-q a c-z-x-v-c\n",
            "**w-r-e-q-w f-a c-z-x-v-c\n",
            "**w-r-e-q-w s c-z-x-v-c\n",
        ]
        with patch("builtins.open", mock_open(read_data="data")):
            result = self.parser.generate_sequences("", add_original=False)

        self.assertEqual(result, sequences)

    def test_generate_sequences_add_mutation_to_original_false(self):
        sequences = [
            "GCATGC\n** a x-z-v-c-x\n",
            "GCATGCAT\n**e f v-c-x-z-v\n",
            "ATGCATGCATGC\n**q-r-e-w-q a x-z-v-c-x\n",
            "TGCATGCATGC\n**r-e-w-q-r f-a x-z-v-c-x\n",
            "ATGCATG\n**q-r-e-w-q s c\n",
            "TGACTG\n** a c-z-x-v-c\n",
            "TGACTGAC\n**r f x-v-c-z-x\n",
            "ACTGACTGACTG\n**q-w-r-e-q a c-z-x-v-c\n",
            "CTGACTGACTG\n**w-r-e-q-w f-a c-z-x-v-c\n",
            "CTGACTGACTG\n**w-r-e-q-w s c-z-x-v-c\n",
        ]
        with patch("builtins.open", mock_open(read_data="data")):
            result = self.parser.generate_sequences("", add_mutation_to_original=False)

        self.assertEqual(result, sequences)

    def test_retrive_sequence(self):
        sequence = "** a x-z-v-c-x"
        tuple_sequence = ("", "a", "xzvcx")

        result = ExtendedParserVcf.retrive_sequence(sequence)

        self.assertEqual(result, tuple_sequence)

    def test_retrive_string_sequence(self):
        sequence = "** a x-z-v-c-x"
        tuple_sequence = "axzvcx"

        result = ExtendedParserVcf.retrive_string_sequence(sequence)

        self.assertEqual(result, tuple_sequence)
