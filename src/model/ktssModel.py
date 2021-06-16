# -*- coding: utf-8 -*-

import functools
import itertools
import json
import logging
import operator
from typing import OrderedDict, Union

from sortedcontainers import SortedDict, SortedSet
from src.argumentParser.abstractArguments import AbstractModelArguments
from src.logging.tqdmLoggingHandler import TqdmLoggingHandler
from src.model.abstractModel import AbstractModel
from src.model.ktssValidation import KTSSValidator
from src.parser.extendedParser import ExtendedParserVcf
from tqdm import tqdm


class KTSSModel(AbstractModel, AbstractModelArguments):
    """Model that creates a DFA from a list of sequences using the ktss method

    Parameters
    ----------
    save_path: bool = False
        The path where the model will be saved
    restore_path: bool = False
        The path from the model ill be loaded
    parser: ParserVcf = ExtendedParserVcf
        Parser that will generate the data
    """

    _arguments: list = [
        {
            "key": "k",
            "name": "k",
            "help": "k value for ktss model",
            "default": 3,
            "type": int,
            "function_argumemnt": {"k_value": "k"},
        },
        {
            "key": "ktss_nas",
            "name": "ktss-not-allowed-segments",
            "help": "Create not allowed segments",
            "function_argumemnt": {"not_allowed_segements": "ktss_nas"},
            "action": "store_true",
        },
    ]
    """ Arguments that will be used by command line """

    _trainer_arguments: dict = {
        "k_value": "k",
        "not_allowed_segements": "get_not_allowed_segements",
    }
    """ Mapping between command line arguments and function arguments of the
    **trainer** method """

    def __init__(
        self,
        save_path=False,
        restore_path=False,
        parser=ExtendedParserVcf,
        tester=KTSSValidator,
    ):
        super().__init__(save_path=save_path, restore_path=restore_path)
        self._model = False
        self.parser_class = parser
        self.tester_class = tester
        self.trainer_name = "ktt"

    def state_in_list(self, state: tuple, lst: list) -> bool:
        """Checks if a state is in a list of states

        Parameters
        ----------
        state: tuple
            State to be searched (s0, s1, s2)
        lst: list
            list where the state will be searched

        Returns
        -------
        True if the state is in the list, if not false
        """
        for state_list in lst:
            if (
                state_list[0] == state[0]
                and state_list[1] == state[1]
                and state_list[2] == state[2]
            ):
                return True
        return False

    def generate_sigma(self, alphabet: Union[list, set], k: int) -> Union[list, set]:
        """Generates sigma for a k given, for example, if k = 2 and alphabet = [A, B]
        it generates [A, B, AA, AB, BA, BB]

        Parameters
        ----------
        alphabet: list, set
            alphabet to generate sigma
        k: int
            length of the maximum strings of sigma

        Returns
        -------
        sigma
        """
        words = alphabet.copy()
        for _ in range(k - 1):
            words.update(
                [pair[0] + pair[1] for pair in itertools.product(words, alphabet)]
            )

        return words

    def get_prefix(self, string: str, k: int) -> str:
        """Returns the prefix of a string

        Parameters
        ----------
        string: str
            Strong where the prefix will be obtained
        k: int
            Length of the prefix - 1

        Returns
        -------
        The prefix of the string
        """
        return string[: k - 1]

    def get_suffix(self, string: str, k: int) -> str:
        """Returns the suffix of a string

        Parameters
        ----------
        string: str
            String where the suffix will be obtained
        k: int
            Length of the suffix - 1

        Returns
        -------
        The suffix of the string
        """
        if k >= len(string):
            return string
        return string[len(string) - k + 1 :]

    def get_infixes(self, string: str, k: int) -> list:
        """Generates all the possible infixes of a string for a lengt k given

        Parameters
        ----------
        string: str
            String where the infixes will be obtained
        k: int
            Length of the infixes

        Rreturns
        --------
        List of all the infixes
        """
        if k >= len(string):
            return [string]
        ini = 0
        fin = k
        result = []
        while fin <= len(string):
            result.append(string[ini:fin])
            ini += 1
            fin += 1
        return result

    def add_transition(
        self,
        transitions: Union[OrderedDict, dict],
        from_state: str,
        symbol: str,
        to_state: str,
    ) -> Union[OrderedDict, dict]:
        """Append a transition into an ordered dict that represent the transitions

        Parameters
        ----------
        transitions: OrderedDict, dict
            Ordered dict that contains transitions
        from_state: str
            String that represent the source state
        symbol: str
            String that represents the transitions symbol
        to_state: str
            String tht represents the destination state

        Returns
        -------
        The updated ordered dictionary with the new transition
        """
        if not transitions.get(from_state):
            transitions[from_state] = SortedDict({})
        transitions[from_state][symbol] = to_state

        return transitions

    def training(
        self,
        samples: Union[list, tuple],
        k: int,
        get_not_allowed_segements: bool = False,
    ) -> Union[OrderedDict, dict]:
        """Generates a ktss model from the samples and a k given

        Parameters
        ----------
        samples: list
            List of samples to training the model
        k: int
            Parameter of the ktss model
        get_not_allowed_segements: bool
            If true returns not allowed segements

        Returns
        -------
            {
                "states": states,
                "alphabet": alphabet,
                "transitions": transitions,
                "initial_state": initial_state,
                "final_states": final_states,
                "not_allowed_segments": not_allowed_segments,
            }
        """
        logger = logging.getLogger()
        tqdm_out = TqdmLoggingHandler(logger, level=logging.INFO)
        logging.info("Training model")
        logging.info("Generating alphabet")
        alphabet = SortedSet(functools.reduce(operator.add, samples))

        greater_or_equal_than_k = []
        lower_than_k = SortedSet()

        logging.info(f"Dviding strings into greater or lower than {k}")
        for sample in samples:
            if len(sample) >= k:
                greater_or_equal_than_k.append(sample)
            else:
                lower_than_k.add(sample)

        logging.info("Generating prefixes and suffixes")
        prefixes = lower_than_k.copy()
        suffixes = lower_than_k.copy()
        infixes = SortedSet()
        for sample in tqdm(greater_or_equal_than_k, file=tqdm_out):
            prefixes.add(self.get_prefix(sample, k))
            suffixes.add(self.get_suffix(sample, k))
            infixes.update(self.get_infixes(sample, k))

        if get_not_allowed_segements:
            logging.info(f"Generating sigma with alphabet {alphabet}")
            not_allowed_segments = self.generate_sigma(alphabet, k) - infixes

        q = [""]
        s = SortedDict({})
        q0 = ""

        logging.info("Generating states from prefixes")
        for prefix in tqdm(prefixes, file=tqdm_out):
            self.add_transition(s, "", prefix[0], prefix[0])
            for char_index in range(len(prefix)):
                q.append([prefix[: char_index + 1]])

                self.add_transition(
                    s, prefix[:char_index], prefix[char_index], prefix[: char_index + 1]
                )

        logging.info("Generating states from infixes")
        for infix in tqdm(infixes, file=tqdm_out):
            q.append([infix[: k - 1], infix[2:k]])
            self.add_transition(s, infix[: k - 1], infix[k - 1], infix[1:k])

        logging.info("Remove repeated and empty states")
        q = SortedSet(
            map(
                lambda x: x if type(x) == str else x[0],
                filter(lambda x: type(x) == str or all(x), q),
            )
        )

        self._model = {
            "states": q,
            "alphabet": alphabet,
            "transitions": s,
            "initial_state": q0,
            "final_states": suffixes,
        }

        if get_not_allowed_segements:
            self._model["not_allowed_segments"] = not_allowed_segments

        logging.info("Training finalized\n")
        return self._model

    @property
    def parser(self):
        return self.parser_class

    @property
    def tester(self):
        return self.tester_class

    @property
    def trainer(self):
        return self.training

    @property
    def model(self):
        return self._model

    def saver(self):
        super().saver()

        with open(f"{self.save_path}ktss-model.json", "w") as outfile:
            model_for_json = {
                "states": list(self.model["states"]),
                "alphabet": list(self.model["alphabet"]),
                "transitions": self.model["transitions"],
                "initial_state": self.model["initial_state"],
                "final_states": list(self.model["final_states"]),
            }

            if self.model.get("not_allowed_segments", False):
                model_for_json["not_allowed_segments"] = list(
                    self.model["not_allowed_segments"]
                )
            json.dump(model_for_json, outfile)

    def loader(self):
        super().loader()

        with open(self.save_path) as json_file:
            model = json.load(json_file)

            model = {
                "states": SortedSet(model["states"]),
                "alphabet": SortedSet(model["alphabet"]),
                "transitions": SortedDict(model["transitions"]),
                "initial_state": model["initial_state"],
                "final_states": SortedSet(model["final_states"]),
            }

            self._model = model

    def get_training_samples(
        self,
        samples: Union[list, tuple],
        has_original: bool = False,
        get_original: bool = False,
    ) -> list:
        """Generates the samples training data from the total of samples

        Parameters
        ----------
        samples: list, tuple
            Samples

        Returns
        -------
        List of samples for training
        """
        return AbstractModel.get_samples(
            samples, self.retrive_string_sequence, has_original, get_original
        )

    def get_test_samples(
        self,
        samples: Union[list, tuple],
        has_original: bool = False,
        get_original: bool = True,
    ) -> list:
        """Generates the samples test data from the total of samples

        Parameters
        ----------
        samples: list, tuple
            Samples

        Returns
        -------
        List of samples for test
        """
        return AbstractModel.get_samples(
            samples, self.retrive_sequence, has_original, get_original
        )

    def get_samples(self, path: str) -> zip:
        with open(path) as samples_file:
            lines = samples_file.readlines()
            original_lines = [
                lines[line] for line in range(len(lines)) if line % 2 == 0
            ]
            parsed_lines = [lines[line] for line in range(len(lines)) if line % 2 == 1]

            zipped_lines = list(zip(original_lines, parsed_lines))
        
        return zipped_lines
