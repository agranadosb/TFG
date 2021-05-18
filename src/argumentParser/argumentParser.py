

from src.constants.constants import (EXTENDED_PARSER_CODE, KTSS_MODEL,
                                     PARSER_MODEL_OPERATION, PARSER_OPERATION)


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
        super().__init__(description)
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
        self.add_argument(
            {
                "key": "p",
                "name": "parser",
                "help": f"Parser to use: extended -> {EXTENDED_PARSER_CODE}",
                "default": EXTENDED_PARSER_CODE,
                "type": str,
                "choices": [EXTENDED_PARSER_CODE],
                "function_argumemnt": {"parser": "parser"},
            }
        )
        self.add_argument(
            {
                "key": "s",
                "name": "save",
                "help": "Folder where save results",
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
