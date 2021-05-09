# -*- coding: utf-8 -*-

import os

from src.constants.constants import (
    EXTENDED_PARSER_CODE,
    KTSS_MODEL,
    LOWER_PARSER_CODE,
    MODEL_OPERATION,
    PARSER_OPERATION,
    SIMPLIFIED_PARSER_CODE,
)
from src.model.ktssModel import KTSSModel
from src.parser.extendedParser import ExtendedParserVcf
from src.parser.lowerParser import LowerParserVcf
from src.parser.simplifiedParser import SimplifiedParserVcf

dir_path = os.path.dirname(os.path.realpath(__file__))


parsers = {
    EXTENDED_PARSER_CODE: ExtendedParserVcf,
    LOWER_PARSER_CODE: LowerParserVcf,
    SIMPLIFIED_PARSER_CODE: SimplifiedParserVcf,
}


def run(
    vcf_path,
    fasta_path,
    model=KTSS_MODEL,
    operation=PARSER_OPERATION,
    parser=SIMPLIFIED_PARSER_CODE,
    result_folder=f"{dir_path}/example/",
):
    if PARSER_OPERATION in operation:
        parser_lower = parsers[parser](vcf_path, fasta_path)

        parser_lower.generate_sequences(result_folder, add_original=False)

    if MODEL_OPERATION in operation:
        pass
