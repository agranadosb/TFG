# -*- coding: utf-8 -*-

import json
import logging
import os

from src.constants.constants import (
    EXTENDED_PARSER_CODE,
    KTSS_MODEL,
    LOWER_PARSER_CODE,
    PARSER_MODEL_OPERATION,
    PARSER_OPERATION,
    SIMPLIFIED_PARSER_CODE,
)
from src.model.ktssModel import KTSSModel
from src.model.ktssValidation import KtssValidator
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

validators = {KTSS_MODEL: KtssValidator}


def run(
    vcf_path,
    fasta_path,
    model_type=KTSS_MODEL,
    operation=PARSER_OPERATION,
    parser=SIMPLIFIED_PARSER_CODE,
    result_folder=f"{dir_path}/example/",
    parser_prefix=20,
    parser_suffix=20,
    test_ratio=0.95,
):
    logging.basicConfig(
        format="%(asctime)s %(levelname)-8s %(message)s",
        level=logging.INFO,
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    if PARSER_MODEL_OPERATION == operation:
        model = models[model_type](save_path=result_folder)

        # Get and parse the data from a file
        parser = model.parser() or parsers[parser]
        parser_engine = parser(vcf_path, fasta_path)
        parser_engine.generate_sequences(
            result_folder,
            add_original=False,
            prefix_length=parser_prefix,
            suffix_length=parser_suffix,
        )

        # Generate the samples to build the model
        with open(f"{result_folder}{parser_engine.default_filename()}") as samples_file:
            lines = samples_file.readlines()
            training_length = int(len(lines) * test_ratio)

            training_samples = model.get_training_samples(lines[:training_length])
            test_samples = model.get_test_samples(lines[training_length:])

        # Train the model
        model.trainer()(training_samples, min(parser_prefix, parser_suffix))
        model.saver()

        # Test the model
        validator = validators[model_type](model.model, parser=model.parser())
        distances = validator.generate_distances(test_samples)

        with open(
            f"{result_folder}{model.trainer_name}-distances.json", "w"
        ) as outfile:
            json.dump(distances, outfile)
        return

    if PARSER_OPERATION in operation:
        parser_method = parsers[parser](vcf_path, fasta_path)
        parser_method.generate_sequences(
            result_folder,
            add_original=False,
            prefix_length=parser_prefix,
            suffix_length=parser_suffix,
        )
