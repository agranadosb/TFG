# -*- coding: utf-8 -*-

import json
import logging
import os
from random import shuffle

from src.argumentParser.argumentParser import ArgumentParser
from src.constants.constants import (EXTENDED_PARSER_CODE, KTSS_MODEL,
                                     MUTATION_PARSER_CODE,
                                     PARSER_MODEL_OPERATION, PARSER_OPERATION)
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
for i in KTSSModel._arguments + ParserVcf._arguments + KTSSValidator._arguments:
    parser.add_argument(i)


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
        self.result_folder = parse_route(result_folder)
        self.model_type = model_type
        self.operation = operation
        parser = parsers[parser]

        self.model = models[model_type](
            save_path=result_folder, parser=parser, tester=validators[model_type]
        )

        self.parser_prefix = parser_prefix
        self.parser_suffix = parser_suffix

        self.save_distances = save_distances
        self.test_ratio = test_ratio
        self.options = kwargs
        self.steps = steps

        self.parser_engine = parser(vcf_path, fasta_path)

    @staticmethod
    def test():
        args = parser.get_function_arguments()

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

        instance = Runner(**parser.get_function_arguments())

        instance.start()

        return instance

    def parse_sequences(self):
        self.parser_engine.generate_sequences(
            self.result_folder,
            prefix_length=self.parser_prefix,
            suffix_length=self.parser_suffix,
            **self.parser_engine.get_generate_sequences_arguments(**self.options),
        )

    def _test_model(self, step):
        test_samples = self.model.get_test_samples(
            self.lines[self.training_length :], has_original=True, get_original=True
        )

        validator = self.model.tester(self.model.model, parser=self.parser_engine)

        distances = validator.generate_distances(
            test_samples,
            prefix_length=self.parser_prefix,
            suffix_length=self.parser_suffix,
            **validator.get_generate_distances_arguments(**self.options),
        )

        ratio_error = 0.0
        for sequence in distances:
            for infix in distances[sequence]:
                ratio_error += distances[sequence][infix]

        distances["error"] = ratio_error

        if self.save_distances:
            with open(
                f"{self.result_folder}{self.model.trainer_name}-distances-{step}.json",
                "w",
            ) as outfile:
                json.dump(distances, outfile)

        return ratio_error

    def _train_model(self):
        training_samples = self.model.get_training_samples(
            self.lines[: self.training_length], has_original=True, get_original=False
        )

        self.model.trainer(
            training_samples, **self.model.get_trainer_arguments(**self.options)
        )
        self.model.saver()

    def _genreate_samples(self):
        self.samples = self.model.get_samples(
            f"{self.result_folder}{self.parser_engine._default_filename}"
        )
        self.training_length = int(len(self.samples) / 2 * self.test_ratio) * 2

    def _shuffle_samples(self):
        shuffle(self.samples)

        self.lines = []
        for line in self.samples:
            self.lines += [line[0], line[1]]

    def train_and_test_model(self):
        total_error = 0.0
        step_ratio = 1 / self.steps
        self._genreate_samples()

        for step in range(self.steps):
            logging.info("###########################################")
            logging.info(f"Validating step {step}")
            logging.info("###########################################")

            self._shuffle_samples()
            self._train_model()
            error_model = self._test_model(step)

            total_error += error_model * step_ratio

        accuracy = (1 - total_error) * 100

        logging.info("###########################################")
        logging.info(f"Model accuracy: {accuracy:.10f} %")
        logging.info("###########################################")

        return accuracy

    def start(self):
        if PARSER_MODEL_OPERATION == self.operation:
            self.parse_sequences()
            accuracy = self.train_and_test_model()
            return accuracy

        if PARSER_OPERATION in self.operation:
            self.parse_sequences()
