# -*- coding: utf-8 -*-

import logging
import os

from src.constants.constants import (
    EXTENDED_PARSER_CODE,
    KTSS_MODEL,
    LOWER_PARSER_CODE,
    MODEL_OPERATION,
    PARSER_MODEL_OPERATION,
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

models = {KTSS_MODEL: KTSSModel}


def run(
    vcf_path,
    fasta_path,
    model=KTSS_MODEL,
    operation=PARSER_OPERATION,
    parser=SIMPLIFIED_PARSER_CODE,
    result_folder=f"{dir_path}/example/",
):
    logging.basicConfig(
        format="%(asctime)s %(levelname)-8s %(message)s",
        level=logging.INFO,
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    if PARSER_MODEL_OPERATION == operation:
        model = models[model](save_path=result_folder)
        parser = model.parser() or parsers[parser]
        parser_engine = parser(vcf_path, fasta_path)
        parser_engine.generate_sequences(result_folder, add_original=False)

        with open(f"{result_folder}{parser_engine.default_filename()}") as samples_file:
            samples = list(
                map(
                    lambda sample: parser_engine.retrive_string_sequence(sample),
                    samples_file,
                )
            )

        model.trainer()(samples, 20)
        model.saver()
        return

    if PARSER_OPERATION in operation:
        parser_method = parsers[parser](vcf_path, fasta_path)
        parser_method.generate_sequences(result_folder, add_original=False)

    if MODEL_OPERATION in operation:
        model = models[model](save_path=result_folder)
        model.trainer()(["abba", "aaabba", "bbaaa", "bba"], 2)
        model.saver()
