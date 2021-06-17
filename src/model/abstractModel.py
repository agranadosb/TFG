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

        self.raw_samples = list(zip(original_lines, parsed_lines))
        self.shuffle_samples()

        self.training_length = int(len(self.samples) / 2 * test_ratio) * 2

        return self.samples

    def shuffle_samples(
        self,
        has_original: bool = True,
    ):
        """Shuffle the samples.

        Parameters
        ----------
        has_original: bool
            Specifies if the file has the original sequence.

        Returns
        -------
        Samples that have been shuffled.
        """
        shuffle(self.raw_samples)
        self.samples = []
        for line in self.raw_samples:
            if has_original:
                self.samples += [line[0], line[1]]
            else:
                self.samples += line

        return self.samples

    def get_training_samples(
        self,
        has_original: bool = True,
        get_original: bool = False,
    ) -> list:
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
        return AbstractModel._get_samples(
            self.samples[: self.training_length],
            self._retrive_string_sample,
            has_original,
            get_original,
        )

    def get_test_samples(
        self,
        has_original: bool = True,
        get_original: bool = True,
    ) -> list:
        """Generates the samples test data from the total of samples.

        Parameters
        ----------
        has_original: bool
            Specifies if the file has the original sequence.
        get_original: bool
            Specifies if the seqeunce that is wanted to be filtered is the original or
            not.

        Returns
        -------
        List of samples for test
        """
        return AbstractModel._get_samples(
            self.samples[self.training_length :],
            self._retrive_sample,
            has_original,
            get_original,
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

    @staticmethod
    def _filter_samples(
        samples: Union[list, tuple],
        has_original: bool = False,
        get_original: bool = False,
    ) -> Union[list, tuple]:
        """If the samples have been parsed with the flag `write_original` to True, this
        method only takes the lines that are in position odd and even, depends if
        `get_orignal` is set True or False.

        This is because the parser creates a file (with the flag `write_original` to
        True) which is an original sequence and parsed sequence per each pair of lines
        in the file.

        It means that the file with the sequences will look like this:

            Original_1
            Parsed_1
            Original_2
            Parsed_2
            ...
            Original_N
            Parsed_N

        So the list will be:

        ```python
            ["Original_1", "Parsed_1", ..., "Original_N", "Parsed_N"]
        ```

        That's the reason for filter depending on the parity of the index of the sample.

        Parameters
        ----------
        samples: list, tuple
            Samples to be filtered.
        has_original: bool
            Specifies if the file has the original sequence.
        get_original: bool
            Specifies if the seqeunce that is wanted to be filtered is the original or
            not.

        Returns
        -------
        List of the filtered samples.
        """
        if has_original:
            index_allowed = 1
            if get_original:
                index_allowed = 0
            samples = [
                samples[i] for i in range(len(samples)) if i % 2 == index_allowed
            ]

        return samples

    @classmethod
    def _get_samples(
        cls,
        samples: Union[list, tuple],
        method: Callable,
        has_original: bool = False,
        get_original: bool = True,
    ) -> list:
        """Gets a list of samples and applies a method to that samples.

        Parameters
        ----------
        samples: list, tuple
            List of samples.
        methods: function
            Function to apply to the samples.
        has_original: bool
            Specifies if the file has the original sequence.
        get_original: bool
            Specifies if the seqeunce that is wanted to be filtered is the original or
            not.

        Returns
        -------
        The list with the samples with the method applied on them.
        """
        samples = cls._filter_samples(samples, has_original, get_original)

        if get_original:
            method = lambda sample: sample

        return list(
            map(
                lambda sample: method(sample.rstrip()),
                samples,
            )
        )
