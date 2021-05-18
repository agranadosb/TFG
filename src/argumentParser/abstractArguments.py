from abc import ABC, abstractmethod
from argparse import ArgumentParser
from typing import Any


class AbstractModelArguments(ABC):
    """Class that allows to implenet command line arguments with a little configuration.

    Each argument is defined as a options in a dict into an arguments list in a model
    class with this required keys:

    - **name**: name of the argument
    - **key**: key of the argument
    - **function_argument**: mapping that maps the name argument with the function
    argument which it will be called

    The other dict keys will be the arguments of argparse add_arguments function.

    It implements functions to get data for methods of models.

    To get the arguments for the trainer method, it's necessary to use the method
    get_trainer_arguments. So, an example could be:

        model = Model(...)
        samples = [...]
        arguments_trainer = {...}
        model.trainer(samples, model.get_trainer_arguments(**arguments_trainer))
    """

    @abstractmethod
    def __init__(self, description: str):
        self.parser: ArgumentParser = ArgumentParser(description=description)
        self.args_map: dict = {}

    @property
    @abstractmethod
    def arguments() -> list:
        return []

    @property
    @abstractmethod
    def trainer_arguments() -> dict:
        return {}

    def get_trainer_arguments(self, **kwargs) -> dict:
        """Generates arguments for the trainer method

        Returns
        -------
        Dictionary with the arguments for the trainer method
        """
        return {
            self.trainer_arguments[key]: value
            for key, value in kwargs.items()
            if key in self.trainer_arguments
        }

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
