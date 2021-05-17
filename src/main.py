# -*- coding: utf-8 -*-

import json
import logging
import os

from src.constants.constants import (
    EXTENDED_PARSER_CODE,
    KTSS_MODEL,
    PARSER_MODEL_OPERATION,
    PARSER_OPERATION,
)
from src.model.ktssModel import KTSSModel
from src.model.ktssValidation import KTSSValidator
from src.parser.extendedParser import ExtendedParserVcf

dir_path = os.path.dirname(os.path.realpath(__file__))


parsers = {EXTENDED_PARSER_CODE: ExtendedParserVcf}

models = {KTSS_MODEL: KTSSModel}

validators = {KTSS_MODEL: KTSSValidator}


def run(
    vcf_path,
    fasta_path,
    model_type=KTSS_MODEL,
    operation=PARSER_OPERATION,
    parser=EXTENDED_PARSER_CODE,
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

        """ TODO: La cadena va a medir con acgt, osea, sin anotar. Hay que ver como
        serán los inputs del validador. Se podría hacer que obtuviera una cadena sin
        anotar y dividiera por 20 (el prefijo que toque) por delante y por detrás para
        generar la cadena para validar"""
        """ TODO: k TIENE QUE SER MENOR A 20, POR EJEMPLO l = 20 y k = 4 (k fijo) """
        # Generate the samples to build the model
        with open(f"{result_folder}{parser_engine.default_filename()}") as samples_file:
            lines = samples_file.readlines()
            training_length = int(len(lines) * test_ratio)

            training_samples = model.get_training_samples(lines[:training_length])
            test_samples = model.get_test_samples(lines[training_length:])

        # Train the model
        model.trainer()(training_samples, 3)
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
