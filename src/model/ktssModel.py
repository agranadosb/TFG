# -*- coding: utf-8 -*-

import functools
import itertools
import json
import logging
import operator
from typing import OrderedDict, Union

from sortedcontainers import SortedDict, SortedList, SortedSet
from src.argumentParser.abstractArguments import AbstractModelArguments
from src.logging.tqdmLoggingHandler import TqdmLoggingHandler
from src.model.abstractModel import AbstractModel
from src.model.ktssValidation import KTSSValidator
from src.model.ktssViterbi import KTSSViterbi
from src.parser.extendedParser import ExtendedParserVcf
from tqdm import tqdm


class KTSSModel(AbstractModel, AbstractModelArguments):
    """Model that creates a DFA from a list of sequences using the ktss method.

    Check [KTSS information](https://riunet.upv.es/bitstream/handle/10251/47562/L%c3%b3pez%20-%20Inference%20of%20k-Testable%20Directed%20Acyclic%20Graph%20Languages.pdf?sequence=1&isAllowed=y)

    Parameters
    ----------
    save_path: bool = False
        The path where the model will be saved.
    restore_path: bool = False
        The path from the model ill be loaded.
    parser: ParserVcf = ExtendedParserVcf
        Parser that will generate the data.
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
    """Arguments that will be used by command line."""

    _trainer_arguments: dict = {
        "k_value": "k",
        "not_allowed_segements": "get_not_allowed_segements",
    }
    """Mapping between command line arguments and function arguments of the
    **trainer** method."""

    def __init__(
        self,
        save_path=False,
        restore_path=False,
        parser=ExtendedParserVcf,
        tester=KTSSViterbi,
    ):
        super().__init__(save_path=save_path, restore_path=restore_path)
        self._model = False
        self.parser_class = parser
        self.tester_class = tester
        self.trainer_name = "ktt"

    def _generate_sigma(self, alphabet: Union[list, set], k: int) -> Union[list, set]:
        """Generates sigma for a k given, for example, if k = 2 and alphabet = [A, B]
        it generates [A, B, AA, AB, BA, BB].

        Parameters
        ----------
        alphabet: list, set
            alphabet to generate sigma.
        k: int
            length of the maximum strings of sigma.

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

    def _get_prefix(self, string: str, k: int) -> str:
        """Returns the prefix of a string.

        Parameters
        ----------
        string: str
            Strong where the prefix will be obtained.
        k: int
            Length of the prefix - 1.

        Returns
        -------
        The prefix of the string.
        """
        return string[: k - 1]

    def _get_suffix(self, string: str, k: int) -> str:
        """Returns the suffix of a string.

        Parameters
        ----------
        string: str
            String where the suffix will be obtained.
        k: int
            Length of the suffix - 1.

        Returns
        -------
        The suffix of the string.
        """
        if k >= len(string):
            return string
        return string[len(string) - k + 1 :]

    def _get_infixes(self, string: str, k: int) -> list:
        """Generates all the possible infixes of a string for a lengt k given.

        Parameters
        ----------
        string: str
            String where the infixes will be obtained.
        k: int
            Length of the infixes.

        Rreturns
        --------
        List of all the infixes.
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

    def _add_transition(
        self,
        transitions: Union[OrderedDict, dict],
        counter: Union[OrderedDict, dict],
        from_state: str,
        symbol: str,
        to_state: str,
    ) -> Union[OrderedDict, dict]:
        """Append a transition into an ordered dict that represent the transitions.

        Parameters
        ----------
        transitions: OrderedDict, dict
            Ordered dict that contains transitions.
        counter: OrderedDict, dict
            Ordered dict that contains the number of times that a transition happens.
        from_state: str
            String that represent the source state.
        symbol: str
            String that represents the transitions symbol.
        to_state: str
            String tht represents the destination state.

        Returns
        -------
        The updated ordered dictionary with the new transition.
        """
        if not transitions.get(from_state):
            transitions[from_state] = SortedDict({})
            counter[from_state] = SortedDict({})

        symbol_counter = 1
        if transitions[from_state].get(symbol, False):
            symbol_counter = counter[from_state][symbol] + 1
        transitions[from_state][symbol] = to_state
        counter[from_state][symbol] = symbol_counter

        return transitions

    @staticmethod
    def _generate_probabilities(
        counter: Union[OrderedDict, dict]
    ) -> Union[OrderedDict, dict]:
        """Creates a dict of probabilities from a counter of transitions. This method
        loops for every transition and generates the probability of each transition per
        symbol dividing bby the posible transitions from a state.

        For example, if our counter is:

        ```python
            {
                "": {"a": 2, "b": 2},
                "a": {"b": 1, "a": 1},
                "b": {"b": 2},
                "aa": {"b": 1, "a": 2},
                "bb": {"a": 4},
                "ba": {"a": 1},
                "ab": {"b": 2},
            }
        ```

        This method will return:

        ```python
            {
                "": {"a": 1 / 2, "b": 1 / 2},
                "a": {"b": 1 / 2, "a": 1 / 2},
                "b": {"b": 1},
                "aa": {"b": 1 / 3, "a": 2 / 3},
                "bb": {"a": 1},
                "ba": {"a": 1},
                "ab": {"b": 1},
            }
        ```

        Parameters
        ----------
        counter: OrderedDict, dict
            Ordered dict that contains the number of times that a transition happens.

        Returns
        -------
        The same dictionary but with probabilities as values.
        """
        result = OrderedDict({})
        for state in counter:
            result[state] = OrderedDict({})
            total = sum(counter[state].values())
            for symbol in counter[state]:
                result[state][symbol] = counter[state][symbol] / total
        return result

    def _generate_sequences(
        self, samples: list, k: int, tqdm_out: TqdmLoggingHandler = None
    ) -> dict:
        """Generates a dictionary with suffixes, ifixes and prefixes for a ktss model.

        Parameters
        ----------
        samples: list
            List of samples where the sequences will be extracted.
        k: int
            K parameters of the ktss model.
        tqdm_out: TqdmLoggingHandler = null
            Tqdm object for loop feedback handling.

        Returns
        -------
        ```python
            {
                "prefixes": prefixes,
                "suffixes": suffixes,
                "infixes": infixes,
            }
        ```
        """
        greater_or_equal_than_k = []
        lower_than_k = SortedList()

        logging.info(f"Dviding strings into greater or lower than {k}")
        for sample in samples:
            if len(sample) >= k:
                greater_or_equal_than_k.append(sample)
            else:
                lower_than_k.add(sample)

        logging.info("Generating prefixes and suffixes")
        prefixes = lower_than_k.copy()
        suffixes = SortedSet(lower_than_k.copy())
        infixes = SortedList()
        for sample in tqdm(greater_or_equal_than_k, file=tqdm_out):
            prefixes.add(self._get_prefix(sample, k))
            suffixes.add(self._get_suffix(sample, k))
            infixes.update(self._get_infixes(sample, k))

        return {"prefixes": prefixes, "suffixes": suffixes, "infixes": infixes}

    def _generate_not_allowed_segments(
        self, infixes: list, alphabet: set, k: int
    ) -> set:
        """Generates not allowed segments for a k given, for example, if k = 2 and
        alphabet = [A, B] and infixes = [AA, BB] it generates [A, B, AB, BA].

        Parameters
        ----------
        infixes: list
            List of allowed infixes.
        alphabet: set
            Alphabet to generate not allowed segments.
        k: int
            Length of the maximum strings of not allowed segments.

        Returns
        -------
        Not allowed segments
        """
        logging.info(f"Generating sigma with alphabet {alphabet}")
        return self._generate_sigma(alphabet, k) - set(infixes)

    def _generate_transitions(
        self,
        initial_state: str,
        prefixes: Union[SortedList, list],
        infixes: Union[SortedList, list],
        k: int,
        tqdm_out: TqdmLoggingHandler = None,
    ):
        """Generates the transitions, the associated probabilities and the tates of a
        ktss model.

        Parameters
        ----------
        initial_state: str
            Initial state
        prefixes: SortedList, list
            List of prefixes of the samples.
        infixes: SortedList, list
            List of infixes of the samples.
        k: int
            K parameters of the ktss model.
        tqdm_out: TqdmLoggingHandler = None
            Tqdm object for loop feedback handling.

        Returns
        -------
        ```python
            {
                "transitions": transitions,
                "states": states,
                "probabilities": probabilities,
            }
        ```
        """
        states = [initial_state]
        transitions = SortedDict({})
        counter = SortedDict({})

        logging.info("Generating states from prefixes")
        for prefix in tqdm(prefixes, file=tqdm_out):
            self._add_transition(
                transitions, counter, initial_state, prefix[0], prefix[0]
            )
            for char_index in range(len(prefix)):
                states.append([prefix[: char_index + 1]])

                self._add_transition(
                    transitions,
                    counter,
                    prefix[:char_index] or initial_state,
                    prefix[char_index],
                    prefix[: char_index + 1],
                )

        logging.info("Generating states from infixes")
        for infix in tqdm(infixes, file=tqdm_out):
            states.append([infix[: k - 1], infix[2:k]])
            self._add_transition(
                transitions, counter, infix[: k - 1], infix[k - 1], infix[1:k]
            )

        logging.info("Remove repeated and empty states")
        states = SortedSet(
            map(
                lambda x: x if type(x) == str else x[0],
                filter(lambda x: type(x) == str or all(x), states),
            )
        )

        probabilities = KTSSModel._generate_probabilities(counter)

        return {
            "transitions": transitions,
            "states": states,
            "probabilities": probabilities,
        }

    def _training(
        self,
        samples: Union[list, tuple],
        k: int,
        get_not_allowed_segements: bool = False,
    ) -> Union[OrderedDict, dict]:
        """Generates a ktss model from the samples and a k given.

        Parameters
        ----------
        samples: list
            List of samples to training the model.
        k: int
            Parameter of the ktss model.
        get_not_allowed_segements: bool
            If true returns not allowed segements.

        Returns
        -------
        ```python
            {
                "states": states,
                "alphabet": alphabet,
                "transitions": transitions,
                "probabilities": probabilities,
                "initial_state": initial_state,
                "final_states": final_states,
                "not_allowed_segments": not_allowed_segments,
            }
        ```
        """
        logger = logging.getLogger()
        tqdm_out = TqdmLoggingHandler(logger, level=logging.INFO)
        logging.info("Training model")

        logging.info("Generating alphabet")
        alphabet = SortedSet(functools.reduce(operator.add, samples))

        sequences = self._generate_sequences(samples, k, tqdm_out)
        infixes = sequences["infixes"]
        prefixes = sequences["prefixes"]
        suffixes = sequences["suffixes"]
        initial_state = "1"

        transitions = self._generate_transitions(
            initial_state, prefixes, infixes, k, tqdm_out
        )

        self._model = {
            "states": transitions["states"],
            "alphabet": alphabet,
            "transitions": transitions["transitions"],
            "initial_state": initial_state,
            "final_states": suffixes,
            "probabilities": transitions["probabilities"],
        }

        if get_not_allowed_segements:
            self._model["not_allowed_segments"] = self._generate_not_allowed_segments(
                infixes, alphabet, k
            )

        logging.info("Training finalized\n")
        return self._model

    @property
    def parser(self):
        return self.parser_class

    def _test(self, parser_engine, filename, save_distances, **kwargs):
        validator = self.tester_class(self.model, parser=parser_engine)

        distances = validator.generate_distances(self.get_test_samples())

        if save_distances:
            with open(filename, "w") as outfile:
                json.dump(distances, outfile)

        return distances["error"]

    @property
    def tester(self):
        return self._test

    def _prepare_training(self, **kwargs):
        self._training(
            self.get_training_samples(), **self.get_trainer_arguments(**kwargs)
        )

    @property
    def trainer(self):
        return self._prepare_training

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
                "probabilities": self.model["probabilities"],
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
