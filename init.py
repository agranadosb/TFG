# -*- coding: utf-8 -*-

import argparse

import src.main as main
from src.constants.constants import (
    EXTENDED_PARSER_CODE,
    KTSS_MODEL,
    LOWER_PARSER_CODE,
    MODEL_OPERATION,
    PARSER_MODEL_OPERATION,
    PARSER_OPERATION,
    SIMPLIFIED_PARSER_CODE,
)

parser = argparse.ArgumentParser(description="Executes a parser, a model or both")

parser.add_argument(
    "-o",
    "--operation",
    help=f"Operation to make: parser -> {PARSER_OPERATION} model -> {MODEL_OPERATION} both -> {PARSER_MODEL_OPERATION}",
    default=PARSER_OPERATION,
    type=str,
    choices=[PARSER_OPERATION, PARSER_MODEL_OPERATION, MODEL_OPERATION],
)

parser.add_argument(
    "-p",
    "--parser",
    help=f"Parser to use: simplified -> {SIMPLIFIED_PARSER_CODE} extended -> {EXTENDED_PARSER_CODE} lower -> {LOWER_PARSER_CODE}",
    default=SIMPLIFIED_PARSER_CODE,
    type=str,
    choices=[SIMPLIFIED_PARSER_CODE, EXTENDED_PARSER_CODE, LOWER_PARSER_CODE],
)

parser.add_argument(
    "-m",
    "--model",
    help=f"Model to use: ktss -> KTSS_MODEL",
    default=KTSS_MODEL,
    type=str,
    choices=[KTSS_MODEL],
)

parser.add_argument(
    "-vcf",
    help="Route to vcf file",
    type=str,
    required=True,
)

parser.add_argument(
    "-fasta",
    help="Route to fasta file",
    type=str,
    required=True,
)

parser.add_argument(
    "-s",
    "--save",
    help="Folder where save results",
    type=str,
)

if __name__ == "__main__":
    args = parser.parse_args()
    main.run(
        args.vcf,
        args.fasta,
        model=args.model,
        operation=args.operation,
        parser=args.parser,
        result_folder=args.save,
    )
