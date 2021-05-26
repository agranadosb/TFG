# -*- coding: utf-8 -*-

import json
import logging
import os
from random import shuffle

from src.argumentParser.argumentParser import ArgumentParser
from src.constants.constants import (
    EXTENDED_PARSER_CODE,
    KTSS_MODEL,
    MUTATION_PARSER_CODE,
    PARSER_MODEL_OPERATION,
    PARSER_OPERATION,
)
from src.model.ktssModel import KTSSModel
from src.model.ktssValidation import KTSSValidator
from src.parser.extendedParser import ExtendedParserVcf
from src.parser.mutationParser import MutationParser
from src.parser.parserVcf import ParserVcf
from src.utils.folders import parse_route

dir_path = os.path.dirname(os.path.realpath(__file__))


parsers = {
    EXTENDED_PARSER_CODE: ExtendedParserVcf,
    MUTATION_PARSER_CODE: MutationParser,
}

models = {KTSS_MODEL: KTSSModel}

validators = {KTSS_MODEL: KTSSValidator}
parser = ArgumentParser("Executes a parser or executes a parser and a model")

""" TODO: Hacer que se detecten automáticamente todos las clases de argumentos """
for i in KTSSModel.arguments + ParserVcf.arguments + KTSSValidator.arguments:
    parser.add_argument(i)


class Runner(object):
    @staticmethod
    def run():
        Runner.start(**parser.get_function_arguments())

    @staticmethod
    def start(
        vcf_path,
        fasta_path,
        model_type=KTSS_MODEL,
        operation=PARSER_OPERATION,
        parser=EXTENDED_PARSER_CODE,
        result_folder=f"{dir_path}/example/",
        parser_prefix=20,
        parser_suffix=20,
        test_ratio=0.95,
        steps=10,
        **kwargs,
    ):
        """Function that runs a model with a parser"""
        logging.basicConfig(
            format="%(asctime)s %(levelname)-8s %(message)s",
            level=logging.INFO,
            datefmt="%Y-%m-%d %H:%M:%S",
        )

        result_folder = parse_route(result_folder)

        if PARSER_MODEL_OPERATION == operation:
            # Get and parse the data from a file
            parser = parsers[parser]

            model = models[model_type](save_path=result_folder, parser=parser)
            parser_engine = parser(vcf_path, fasta_path)
            parser_engine.generate_sequences(
                result_folder,
                prefix_length=parser_prefix,
                suffix_length=parser_suffix,
                **parser_engine.get_generate_sequences_arguments(**kwargs),
            )

            # Generate the samples to build the model
            with open(
                f"{result_folder}{parser_engine.default_filename}"
            ) as samples_file:
                lines = samples_file.readlines()
                original_lines = [
                    lines[line] for line in range(len(lines)) if line % 2 == 0
                ]
                parsed_lines = [
                    lines[line] for line in range(len(lines)) if line % 2 == 1
                ]

                zipped_lines = list(zip(original_lines, parsed_lines))

            for step in range(steps):
                logging.info("###########################################")
                logging.info(f"Validating step {step}")
                logging.info("###########################################")

                shuffle(zipped_lines)
                lines = []
                for line in zipped_lines:
                    lines += [line[0], line[1]]

                training_length = int(len(lines) / 2 * test_ratio) * 2

                training_samples = model.get_training_samples(
                    lines[:training_length], has_original=True, get_original=False
                )
                test_samples = model.get_test_samples(
                    lines[training_length:], has_original=True, get_original=True
                )

                # Train the model
                model.trainer(training_samples, **model.get_trainer_arguments(**kwargs))
                model.saver()

                # Test the model
                validator = validators[model_type](model.model, parser=parser)

                distances = validator.generate_distances(
                    test_samples,
                    prefix_length=parser_prefix,
                    suffix_length=parser_suffix,
                    **validator.get_generate_distances_arguments(**kwargs),
                )

                with open(
                    f"{result_folder}{model.trainer_name}-distances-{step}.json", "w"
                ) as outfile:
                    json.dump(distances, outfile)
            return

        if PARSER_OPERATION in operation:
            parser_method = parsers[parser](vcf_path, fasta_path)
            parser_method.generate_sequences(
                result_folder,
                add_original=True,
                prefix_length=parser_prefix,
                suffix_length=parser_suffix,
            )
