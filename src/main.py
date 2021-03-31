# -*- coding: utf-8 -*-

import os

from src.vcf.vcfReader import VcfMutationsReader
from src.model.parserVcf import ParserVcf
from src.tests.vcfReaderTest import test_get_seq_by_chr_pos


dir_path = os.path.dirname(os.path.realpath(__file__))


def run(test = False):
    if test:
        test_get_seq_by_chr_pos()
        return
    reader = ParserVcf(
        f'{dir_path}/example/datosR1.vcf',
        f'{dir_path}/example/hg19.fa.gz'
    )

    index = 0
    for i in reader.get_vcf():
        if index == 13:
            break
        print(i.REF)
        print(reader.get_vcf_reader().get_sequence(i.CHROM, i.REF, i.POS, 20, 20))
        index += 1
