# -*- coding: utf-8 -*-

from unittest import TestCase

from src.vcf.vcfReader import VcfMutationsReader


class TestVcfMutationsReader(TestCase):
    def create_mockup_vcf(self):
        return VcfMutationsReader(
            '/opt/UPV/TFG/src/tests/vcfTest.vcf',
            '/opt/UPV/TFG/src/tests/test.fa.gz'
        )

    def test_set_fasta_line_length(self):
        reader = self.create_mockup_vcf()

        reader.set_fasta_line_length()

        assert reader.fasta_file_line_length == 50

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
