# -*- coding: utf-8 -*-


from abc import ABC, abstractmethod
from random import shuffle
from typing import Callable, Union

from src.utils.folders import parse_route


class AbstractModel(ABC):
    """This class contains default functionality to generate a model that will have
    functions to train a model and functions to get the generated model.

    This class will need data to train the model. This data will be typically a list of
    pairs (sample, result) or a list of samples that will be used to train a specified
    model.

    The model will be a defined data structure and will be created using a trainer.

    A trainer is a function that gets samples as input (typically a list of pairs or a
    list of samples) and returns the trained model.

    The steps will be:

     - DATA -----> PARSER -----> TRAINER -----> MODEL -----> TESTER [-----> SAVER]
     - DATA -----> PARSER [-----> LOADER] -----> MODEL -----> TESTER -----> RESULT

    So, to implement this class, will be necessary to define the next methods and
    attributes:

     - **parser**: Attribute which is the class of the parser used by the model.
     - **trainer**: Function to train the model.
     - **model**: Attribute which returns the model after training.
     - **tester**: Function that test the model.
     - **saver**: Method that saves the model in a file.
     - **loader**: Method that loads the model from a file.

    This class also has methods to manage the samples for test and training methods:

     - **get_samples**: Gets samples from a file that has a pair sample, one item per
     line.
     - **shuffle_samples**: Shuffle the samples.
     - **get_training_samples**: Generates the samples training data from the total of
     samples.
     - **get_test_samples**: Generates the samples test data from the total of samples.

    Parameters
    ----------
    save_path: str
        The path where the results will be stored.
    restore_path: str
        Path from where the results will be retrieved.
    """

    @abstractmethod
    def __init__(self, save_path: str, restore_path: str):
        self.save_path = parse_route(save_path)
        self.restore_path = parse_route(restore_path)

    @property
    @abstractmethod
    def parser(self):
        """Attriubte which is the class of the parser used by the model."""
        pass

    @abstractmethod
    def trainer(self):
        """Function to train the model."""
        pass

    @property
    @abstractmethod
    def model(self):
        """Attribute which returns the model after training."""
        pass

    @abstractmethod
    def tester(self):
        """Function that test the model."""
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

    def get_samples(self, path: str, test_ratio: int, has_original: bool = True) -> zip:
        """Gets samples from a file that has a pair sample, one item per line.

        Parameters
        ----------
        path: str
            Path of the file with the samples.
        test_ratio: int
            Ratio of training and test samples
        has_original: bool
            Specifies if the file has the original sequence.

        Returns
        -------
        Samples in a list of pairs.
        """
        with open(path) as samples_file:
            lines = samples_file.readlines()

        if not has_original:
            return lines

        original_lines = [lines[line] for line in range(len(lines)) if line % 2 == 0]
        parsed_lines = [lines[line] for line in range(len(lines)) if line % 2 == 1]

        self.samples = list(zip(original_lines, parsed_lines))
        self.training_length = int(len(self.samples) * test_ratio)

        return self.samples

    def shuffle_samples(self):
        """Shuffle the samples.

        Returns
        -------
        Samples that have been shuffled.
        """
        shuffle(self.samples)

    def get_training_samples(self, get_original: bool = True) -> list:
        """Generates the samples training data from the total of samples.

        Parameters
        ----------
        has_original: bool
            Specifies if the file has the original sequence.
        get_original: bool
            Specifies if the seqeunce that is wanted to be filtered is the original or
            not.

        Returns
        -------
        List of samples for training
        """
        samples = self.samples
        if get_original:
            samples = list(
                map(lambda sample: sample[1], self.samples[: self.training_length])
            )

        return AbstractModel._get_samples(
            samples[: self.training_length], self._retrive_string_sample
        )

    def get_test_samples(self) -> list:
        """Generates the samples test data from the total of samples.

        Parameters
        ----------
        has_original: bool
            Specifies if the file has the original sequence.

        Returns
        -------
        List of samples for test
        """
        return AbstractModel._get_samples(
            self.samples[self.training_length :], self._retrive_string_sample
        )

    @property
    def _retrive_string_sample(self) -> Callable:
        """Gets a sequence in a string format and returns it in a different string
        format.

        Parameters
        ----------
        sequence: str
            Sequence in string format.

        Returns
        -------
        The sequence in a string format.
        """
        return self.parser.retrive_string_sequence

    @property
    def _retrive_sample(self) -> Callable:
        """Gets a sequence in a string format and returns it in a tuple format.

        Parameters
        ----------
        sequence: str
            Sequence in string format.

        Returns
        -------
        Sequence in a list format.
        """
        return self.parser.retrive_sequence

    @classmethod
    def _get_samples(cls, samples: Union[list, tuple], method: Callable) -> list:
        """Gets a list of samples and applies a method to that samples.

        Parameters
        ----------
        samples: list, tuple
            List of samples.
        method: function
            Function to apply to the samples.

        Returns
        -------
        The list with the samples with the method applied on them.
        """
        result = []
        for sample in samples:
            if isinstance(sample, (list, tuple)):
                sample = (sample[0].rstrip(), method(sample[1].rstrip()))
            else:
                sample = method(sample.rstrip())
            result.append(sample)

        return result
