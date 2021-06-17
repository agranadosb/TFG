# -*- coding: utf-8 -*-

import logging
import os

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


_parsers = {
    EXTENDED_PARSER_CODE: ExtendedParserVcf,
    MUTATION_PARSER_CODE: MutationParser,
}

_models = {KTSS_MODEL: KTSSModel}

_validators = {KTSS_MODEL: KTSSValidator}
_argument_parser = ArgumentParser("Executes a parser or executes a parser and a model")

""" TODO: Hacer que se detecten automáticamente todos las clases de argumentos """
for i in KTSSModel._arguments + ParserVcf._arguments + KTSSValidator._arguments:
    _argument_parser.add_argument(i)


class Runner(object):
    def __init__(
        self,
        vcf_path,
        fasta_path,
        model_type=KTSS_MODEL,
        operation=PARSER_OPERATION,
        parser=EXTENDED_PARSER_CODE,
        result_folder=f"{dir_path}/example/",
        parser_prefix=20,
        parser_suffix=20,
        test_ratio=0.95,
        save_distances=False,
        steps=10,
        **kwargs,
    ):
        self._result_folder = parse_route(result_folder)
        self._operation = operation
        parser = _parsers[parser]

        self._model = _models[model_type](
            save_path=result_folder, parser=parser, tester=_validators[model_type]
        )

        self._parser_prefix = parser_prefix
        self._parser_suffix = parser_suffix

        self._save_distances = save_distances
        self._test_ratio = test_ratio
        self._options = _argument_parser.get_function_arguments()
        self._steps = steps

        self._parser_engine = parser(vcf_path, fasta_path)

    def parse_sequences(self):
        self._parser_engine.generate_sequences(**self._options)

    def train_and_test_model(self):
        total_error = 0.0
        step_ratio = 1 / self._steps
        self._model.get_samples(
            f"{self._result_folder}{self._parser_engine._default_filename}",
            self._test_ratio,
        )

        for step in range(self._steps):
            logging.info("###########################################")
            logging.info(f"Validating step {step}")
            logging.info("###########################################")

            self._model.shuffle_samples()
            self._model.trainer(**self._options)
            self._model.saver()

            filename = (
                f"{self._result_folder}{self._model.trainer_name}-distances-{step}.json"
            )

            error_model = self._model.tester(
                self._parser_engine,
                self._save_distances,
                filename,
                **self._options,
            )

            total_error += error_model * step_ratio

        accuracy = (1 - total_error) * 100

        logging.info("###########################################")
        logging.info(f"Model accuracy: {accuracy:.10f} %")
        logging.info("###########################################")

        return accuracy

    @staticmethod
    def test():
        args = _argument_parser.get_function_arguments()

        with open("results.txt", "w") as fr:
            print(f"LENGTH\tK\tRATIO\tACCURACY", file=fr)
            for length in range(5, 51, 10):
                for k in range(2, 10):
                    for ratio in range(5, 10):
                        args["test_ratio"] = ratio / 10
                        args["parser_prefix"] = length
                        args["parser_suffix"] = length
                        args["k_value"] = k

                        print(
                            f"{length}\t{k}\t{ratio}\t{Runner.run():.2f}",
                            file=fr,
                        )

    @staticmethod
    def run(logger=True):
        logging.basicConfig(
            format="%(asctime)s %(levelname)-8s %(message)s",
            level=logging.INFO if logger else logging.WARNING,
            datefmt="%Y-%m-%d %H:%M:%S",
        )

        instance = Runner(**_argument_parser.get_function_arguments())

        instance._start()

        return instance

    def _start(self):
        if PARSER_MODEL_OPERATION == self._operation:
            self.parse_sequences()
            return self.train_and_test_model()

        if PARSER_OPERATION in self._operation:
            self.parse_sequences()
