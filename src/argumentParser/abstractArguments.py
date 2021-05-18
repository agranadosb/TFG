from abc import ABC, abstractmethod


class AbstractModelArguments(ABC):
    """Class that allows to implenet command line arguments with a little configuration.

    Each argument is defined as a dict into an arguments list in a model class with this
    required keys:
        - name: name of the argument
        - key: key of the argument
        - function_argument: mapping that maps the name argument with the function
        argument which it will be called

    The other dict keys will be the arguments of argparse add_arguments function.

    It implements functions to get data for methods of models.

    To get the arguments for the trainer method, it's necessary to user ed method
    get_trainer_arguments. So, an example could be:

        model = Model(...)
        samples = [...]
        arguments_trainer = {...}
        model.trainer(samples, model.get_trainer_arguments(**arguments_trainer))
    """

    @property
    @abstractmethod
    def arguments():
        return []

    @property
    @abstractmethod
    def trainer_arguments():
        return {}

    def get_trainer_arguments(self, **kwargs):
        return {
            self.trainer_arguments[key]: value
            for key, value in kwargs.items()
            if key in self.trainer_arguments
        }
