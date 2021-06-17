from argparse import ArgumentParser as DefaultParser
from typing import Any

from src.constants.constants import (
    EXTENDED_PARSER_CODE,
    KTSS_MODEL,
    MUTATION_PARSER_CODE,
    PARSER_MODEL_OPERATION,
    PARSER_OPERATION,
)


class ArgumentParser(object):
    """Class that add arguments to be parsed and setted into run function.

    An argument can be added using the method add_argument, which gets a dict with
    requierd keys:

    - name: name of the argument
    - key: key of the argument
    - function_argument: mapping that maps the name argument with the function
    argument which it will be called

    The other dict keys will be the arguments of argparse add_arguments function.
    """

    def __init__(self, description: str):
        self.parser: DefaultParser = DefaultParser(description=description)
        self.args_map: dict = {}

        """ TODO: Hacer un atributo para añadir los modelos """
        self.add_argument(
            {
                "key": "m",
                "name": "model",
                "help": f"Model to use: ktss -> {KTSS_MODEL}",
                "default": KTSS_MODEL,
                "type": str,
                "choices": [KTSS_MODEL],
                "function_argumemnt": {"model_type": "model"},
            }
        )
        self.add_argument(
            {
                "key": "o",
                "name": "operation",
                "help": f"Operation to make: parser -> {PARSER_OPERATION}, both -> {PARSER_MODEL_OPERATION}",
                "default": PARSER_OPERATION,
                "type": str,
                "choices": [PARSER_OPERATION, PARSER_MODEL_OPERATION],
                "function_argumemnt": {"operation": "operation"},
            }
        )
        """ TODO: Hacer un atributo para añadir los parsers """
        self.add_argument(
            {
                "key": "p",
                "name": "parser",
                "help": f"Parser to use: extended -> {EXTENDED_PARSER_CODE}, mutation type -> {MUTATION_PARSER_CODE}",
                "default": EXTENDED_PARSER_CODE,
                "type": str,
                "choices": [EXTENDED_PARSER_CODE, MUTATION_PARSER_CODE],
                "function_argumemnt": {"parser": "parser"},
            }
        )
        self.add_argument(
            {
                "key": "s",
                "name": "save",
                "help": "Folder where save results",
                "required": True,
                "type": str,
                "function_argumemnt": {"result_folder": "save"},
            }
        )
        self.add_argument(
            {
                "key": "p_p",
                "name": "parser_prefix",
                "help": "Length of the sequence prefix",
                "default": 20,
                "type": int,
                "function_argumemnt": {"parser_prefix": "parser_prefix"},
            }
        )
        self.add_argument(
            {
                "key": "p_s",
                "name": "parser_suffix",
                "help": "Length of the sequence suffix",
                "default": 20,
                "type": int,
                "function_argumemnt": {"parser_suffix": "parser_suffix"},
            }
        )
        self.add_argument(
            {
                "key": "r",
                "name": "ratio",
                "help": "Ratio of training and test samples",
                "default": 0.9,
                "type": float,
                "function_argumemnt": {"test_ratio": "ratio"},
            }
        )
        self.add_argument(
            {
                "key": "vcf",
                "name": "vcf",
                "help": "Route to vcf file",
                "required": True,
                "type": str,
                "function_argumemnt": {"vcf_path": "vcf"},
            }
        )
        self.add_argument(
            {
                "key": "fasta",
                "name": "fasta",
                "help": "Route to fasta file",
                "required": True,
                "type": str,
                "function_argumemnt": {"fasta_path": "fasta"},
            }
        )
        self.add_argument(
            {
                "key": "steps",
                "name": "steps",
                "help": "Rounds to execute the validator",
                "default": 10,
                "type": int,
                "function_argumemnt": {"steps": "steps"},
            }
        )
        self.add_argument(
            {
                "key": "sd",
                "name": "save-distances",
                "help": "Save he distances into files per step",
                "function_argumemnt": {"sd": "save_distances"},
                "action": "store_true",
            },
        )

    def add_argument(self, options: dict):
        """Adds an argument to the argument parser

        Parameters
        ----------
        options: dict
            Dictionary that contains the options for the argument
        """
        function_argument = list(options.pop("function_argumemnt").items())[0]
        name = options.pop("name")
        key = options.pop("key")

        self.args_map[function_argument[0]] = function_argument[1]
        self.parser.add_argument(
            f"-{key}",
            f"--{name}",
            **options,
        )

    def get_argument(self, key: str) -> Any:
        """Gets an argument from argument parser

        Parameters
        ----------
        key: str
            Key of the argument

        Returns
        -------
        The value of the argument
        """
        args = self.parse_args()

        value = self.args_map[key]

        return getattr(args, value, False)

    def get_function_arguments(self) -> dict:
        """Creates a dictionary mapping to the function arguments with all the arguments
        of the argument parser

        Returns
        -------
        Mapping of the arguments and the function arguments
        """
        return {key: self.get_argument(key) for key in self.args_map}

    def parse_args(self) -> dict:
        """Parse the command line args

        Returns
        -------
        Dictionary with the arguments
        """
        return self.parser.parse_args()
