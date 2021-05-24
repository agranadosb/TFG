from abc import ABC, abstractmethod

from src.argumentParser.argumentParser import ArgumentParser


class AbstractValidationArguments(ABC, ArgumentParser):
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

    @property
    @abstractmethod
    def arguments() -> list:
        return []

    @property
    @abstractmethod
    def generate_distances_arguments() -> dict:
        return {}

    def get_generate_distances_arguments(self, **kwargs) -> dict:
        """Generates arguments for the generate_distances method

        Returns
        -------
        Dictionary with the arguments for the generate_distances method
        """
        return {
            self.generate_distances_arguments[key]: value
            for key, value in kwargs.items()
            if key in self.generate_distances_arguments
        }


class AbstractParserArguments(ABC, ArgumentParser):
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

    @property
    @abstractmethod
    def arguments() -> list:
        return []

    @property
    @abstractmethod
    def generate_sequences_arguments() -> dict:
        return {}

    def get_generate_sequences_arguments(self, **kwargs) -> dict:
        """Generates arguments for the generate_sequences method

        Returns
        -------
        Dictionary with the arguments for the generate_sequences method
        """
        return {
            self.generate_sequences_arguments[key]: value
            for key, value in kwargs.items()
            if key in self.generate_sequences_arguments
        }


class AbstractModelArguments(ABC, ArgumentParser):
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
