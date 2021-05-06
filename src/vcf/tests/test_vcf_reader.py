# -*- coding: utf-8 -*-

import pathlib
from unittest import TestCase

from src.vcf.vcfReader import VcfMutationsReader


class TestVcfMutationsReader(TestCase):
    def setUp(self) -> None:
        self.reader = VcfMutationsReader(
            "/opt/UPV/TFG/src/vcf/tests/vcfTest.vcf",
            "/opt/UPV/TFG/src/vcf/tests/test.fa.gz",
        )

        return super().setUp()

    def test_set_fasta_line_length(self):
        line_length = 50

        self.assertEqual(self.reader.fasta_file_line_length, line_length)

    def test_set_fasta_file(self):
        path = f"{str(pathlib.Path(__file__).parent.absolute())}/test.fa"

        self.assertEqual(self.reader.fasta_file.name, path)

    def test_set_chromosme_information_first(self):
        chromosme_label = "chr1"
        label_length = len("chr1")
        line_starts = 0
        index_starts = 0

        self.assertTrue(
            self.reader.fatsa_chromosme_information.get(chromosme_label, True)
        )
        self.assertEqual(
            self.reader.fatsa_chromosme_information[chromosme_label][label_length],
            label_length,
        )
        self.assertEqual(
            self.reader.fatsa_chromosme_information[chromosme_label][line_starts],
            line_starts,
        )
        self.assertEqual(
            self.reader.fatsa_chromosme_information[chromosme_label][index_starts],
            index_starts,
        )


""" # -*- coding: utf-8 -*-

from src.vcf.vcfReader import VcfMutationsReader

results = [
    ('', 'T', 'ATCAA'),
    ('', 'TA', 'TCAAT'),
    ('T', 'A', 'TCAAT'),
    ('TA', 'T', 'CAATG'),
    ('TAT', 'C', 'AATGC'),
    ('TATC', 'A', 'ATGCC'),
    ('TATCA', 'A', 'TGCCT'),
    ('GAGGG', 'C', 'ATAAC'),
    ('GAGGG', 'CA', 'TAACA'),
    ('AGGGC', 'A', 'TAACA'),
    ('GGGCA', 'T', 'AACAT'),
    ('GGCAT', 'A', 'ACATC'),
    ('AGGTA', 'T', 'CTAAT'),
    ('GGTAT', 'C', 'TAAT'),
    ('GTATC', 'T', 'AAT'),
    ('TATCT', 'A', 'AT'),
    ('ATCTA', 'A', 'T'),
    ('ATCTA', 'AT', ''),
    ('TCTAA', 'T', '')
]


def test_get_seq_by_chr_pos():
    reader = VcfMutationsReader(
        '/opt/UPV/TFG/src/tests/vcfTest.vcf',
        '/opt/UPV/TFG/src/tests/test.fa.gz'
    )

    for test_vcf_line in zip(reader.get_vcf(), results * 3):
        vcf_line = test_vcf_line[0]
        result = test_vcf_line[1]
        reference_seq_fasta = reader.get_sequence(
            vcf_line.CHROM,
            vcf_line.REF,
            vcf_line.POS,
            5,
            5
        )

        assert result[0] == reference_seq_fasta[0]
        assert vcf_line.REF == reference_seq_fasta[1]
        assert result[2] == reference_seq_fasta[2] """

"""from src.model.parserVcf import ParserVcf


def test_generate_simplified_sequences():
    parser_simplified = ParserVcf(
        '/opt/UPV/TFG/src/tests/vcfTest.vcf',
        '/opt/UPV/TFG/src/tests/test.fa.gz'
    )

    parser_lower = ParserVcf(
        '/opt/UPV/TFG/src/tests/vcfTest.vcf',
        '/opt/UPV/TFG/src/tests/test.fa.gz'
    )

    # TODO: Hay que testear get_simplified_sequence

    parser_simplified.generate_simplified_sequences(
        '/opt/UPV/TFG/src/tests'
    )

    parser_lower.generate_lower_sequences(
        '/opt/UPV/TFG/src/tests'
    )
"""
