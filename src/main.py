# -*- coding: utf-8 -*-

import os

from src.model.parserVcf import ParserVcf

dir_path = os.path.dirname(os.path.realpath(__file__))


def run(test=False):
    parser_lower = ParserVcf(
        f"{dir_path}/example/datosR1.vcf", f"{dir_path}/example/hg19.fa.gz"
    )

    parser_lower.generate_extended_sequences(f"{dir_path}/example/", add_original=False)

    """ index = 0
    for i in reader.get_vcf():
        if index == 13:
            break
        print(i.REF)
        print(reader.get_vcf_reader().get_sequence(i.CHROM, i.REF, i.POS, 20, 20))
        index += 1 """
