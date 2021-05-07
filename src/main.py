# -*- coding: utf-8 -*-

import os

from src.parser.extendedParser import ExtendedParserVcf

dir_path = os.path.dirname(os.path.realpath(__file__))


def run():
    parser_lower = ExtendedParserVcf(
        f"{dir_path}/example/datosR1.vcf", f"{dir_path}/example/hg19.fa.gz"
    )

    parser_lower.generate_sequences(f"{dir_path}/example/", add_original=False)
