# -*- coding: utf-8 -*-


from abc import ABC, abstractmethod
from typing import Callable, Union

from src.utils.folders import parse_route


class AbstractModel(ABC):
    """
    This class contains default functionallity to generate a model that will have
    functions to train and functions to get the generated model.

    This class will need data to train the model. This data will be tipically a list of
    pairs (sample, result) or a list of samples that will be used to train an specified
    model.

    The model will be a defined data structure, and will be created using a trainer.

    A trainer is a function that gets samples as input (tipically a list of pairs or a
    list of samples) and returns the trained model.

    The steps will be:

        DATA -----> PARSER -----> TRAINER -----> MODEL [-----> SAVER]
        DATA -----> PARSER [-----> LOADER] -----> MODEL -----> RESULT

    So, to implement this class, will be necessary to define the next methods and
    attributes:

     - **parser**: Attriubte which is the class of the parser used by the model
     - **trainer**: Attriubte which returns the function to train the model
     - **model**: Attribute which returns the model after training
     - **saver**: Method that save the model in a file
     - **loader**: Method that loads the model from a file
    """

    @abstractmethod
    def __init__(self, save_path: str = None, restore_path: str = None):
        self.save_path = parse_route(save_path)
        self.restore_path = parse_route(restore_path)

    @property
    @abstractmethod
    def parser(self):
        """Attriubte which is the class of the parser used by the model."""
        pass

    @property
    @abstractmethod
    def trainer(self):
        """Attriubte which returns the function to train the model."""
        pass

    @property
    @abstractmethod
    def model(self):
        """Attribute which returns the model after training."""
        pass

    @abstractmethod
    def saver(self):
        """Method that save the model in a file."""
        if not self.save_path:
            raise AttributeError("Saver path is not defined")

        if not self.model:
            raise ValueError("The model is not trained")

    @abstractmethod
    def loader(self):
        """Method that loads the model from a file."""
        if not self.restore_path:
            raise AttributeError("Loader path is not defined")

    @property
    def retrive_string_sequence(self) -> Callable:
        """Gets a string sequence and returns the sequence in a string type.

        Parameters
        ----------
        sequence: str
            Sequence in string format

        Returns
        -------
        The sequence in a string format
        """
        return self.parser.retrive_string_sequence

    @property
    def retrive_sequence(self) -> Callable:
        """Gets a string sequence and returns the sequence in a tuple type.

        Parameters
        ----------
        sequence: str
            Sequence in string format

        Returns
        -------
        Sequence in a list format
        """
        return self.parser.retrive_sequence

    def get_samples(self, samples: Union[list, tuple], method: Callable) -> list:
        """Gets a list of samples and applies a method to that samples.

        Parameters
        ----------
        samples: list, tuple
            List of samples
        methods: function
            Function to apply to the samples

        Returns
        -------
        The list with the samples with the method applied on them
        """
        return list(
            map(
                lambda sample: method(sample.rstrip()),
                samples,
            )
        )
