# -*- coding: utf-8 -*-

from src.vcf.vcfReader import VcfMutationsReader


class TestVcfMutationsReader:
    def create_mockup_vcf(self):
        return VcfMutationsReader(
            '/opt/UPV/TFG/src/tests/vcfTest.vcf',
            '/opt/UPV/TFG/src/tests/test.fa.gz'
        )

    def test_set_fasta_line_length(self):
        reader = self.create_mockup_vcf()

        reader.set_fasta_line_length()

        assert reader.fasta_file_line_length == 50
