import logging
import random
from typing import Callable, Union

from sortedcontainers import SortedDict
from src.argumentParser.abstractArguments import AbstractValidationArguments
from src.dataStructures.dfaStochastic import DFAStochastic
from src.logging.tqdmLoggingHandler import TqdmLoggingHandler
from src.parser.extendedParser import ExtendedParserVcf
from src.parser.parserVcf import ParserVcf
from textdistance import levenshtein
from tqdm import tqdm


class KTSSValidator(AbstractValidationArguments):
    """Operates from a DFA generated from a KTSS model and allows to generate distances
    between original sequences and sequences derived from the original sequences.

    Parameters
    ----------
    model: dict
        Dictionary that contains de DFA model
    parser: ParserVCF = ExtendedParserVcf
        Parser of the model
    """

    _arguments: list = [
        {
            "key": "sep",
            "name": "separator",
            "help": "Specifies the separator between characters of each sequence of the valdiator",
            "default": "-",
            "type": str,
            "function_argumemnt": {"sep": "separator"},
        },
        {
            "key": "min",
            "name": "minimum",
            "help": "If true only returns the minimum value and infix of the all the distances of the valdiator",
            "function_argumemnt": {"min": "minimum"},
            "action": "store_true",
        },
        {
            "key": "aoval",
            "name": "add-original-validator",
            "help": "If true returns the original sequence, if not returns the anotated sequence of the valdiator",
            "function_argumemnt": {"aoval": "add_original_validator"},
            "action": "store_true",
        },
        {
            "key": "amval",
            "name": "add-mutation-validator",
            "help": "If true returns the mutation, if not returns the reference sequence",
            "function_argumemnt": {"amval": "add_mutation_validator"},
            "action": "store_true",
        },
    ]
    """ Arguments that will be used by command line """

    _generate_distances_arguments: dict = {
        "sep": "separator",
        "min": "minimum",
        "aoval": "add_original",
        "amval": "add_mutation",
        "parser_prefix": "prefix_length",
        "parser_suffix": "suffix_length",
    }
    """ Mapping between command line arguments and function arguments of the
    **generate_distances** method """

    def __init__(
        self, model: Union[SortedDict, dict], parser: ParserVcf = ExtendedParserVcf
    ):
        self.model = model
        self.infix_symbols = parser.mutations_symbols

        self._set_mappings(parser)
        self.parser = parser
        self.dfa = DFAStochastic(
            model["states"],
            model["alphabet"],
            model["transitions"],
            model["initial_state"],
            model["final_states"],
            model["probabilities"],
        )

    def generate_distances(self, sequences: Union[list, tuple]) -> SortedDict:
        """Generates all the distances of an infix of a given list of sequences of all
        the possible infixes.

        Parameters
        ---------
        sequences: list
            List of sequences.

        Returns
        -------
        Dictionary that contains the sequence and its distances.
        """
        logger = logging.getLogger()
        tqdm_out = TqdmLoggingHandler(logger, level=logging.INFO)

        result = SortedDict()
        logger.info("Generating validation data")
        total_errors = 0
        total_chars = 0
        for sequence_raw in tqdm(sequences, file=tqdm_out):
            reference_sequence = sequence_raw[0]
            annotated_sequence = sequence_raw[1]

            result[reference_sequence] = self._string_distances(
                self.annotate_sequence(reference_sequence), annotated_sequence
            )

            total_errors += result[reference_sequence]
            total_chars += len(annotated_sequence)

        if total_chars != 0:
            result["error"] = total_errors / total_chars

        return result

    def annotate_sequence(self, sequence: str, separator: str = "") -> str:
        """Gets a string sequence, and annotates it.

        The process consists of getting all the associated symbols from the parser
        mappings (prefix, mutation, and suffix) that can be associated with a symbol
        from the sequence.

        Then, the method filters them by the condition that a transition can be done
        with that symbols from a state of the model DFA.

        If more than one symbol can be associated to a sequence symbol, this methods get
        in a aleatory way a symbol from all the possible symbols.

        For instance, if our sequence is:

        ```python
            "AAA"
        ```

        Our transitions are:

        ```python
            {
                "1": {"a": "2", "b": "3"},
                "2": {"a": "4", "b": "5", "d": "6", "l": "15"},
                "3": {"a": "7", "b": "8", "d": "9"},
                "4": {"z": "10", "x": "11"},
                "6": {"z": "12"},
                "9": {"z": "13", "x": "14"},
                "15": {"z": "16"},
            }
        ```

        And our mappings are:

        ```python
            prefix_map = {"A": "a", "C": "c", "B": "b"}
            mutations_map = {"A": "l", "C": "d"}
            suffix_map = {"A": "z", "C": "x", "B": "c"}
        ```

        We can get two results (with 50% of chance per result):

         - **aaz**
         - **alz**
        """
        result = []
        current_state = self.dfa.initial_state
        for symbol in sequence:
            possible_symbols = []
            mutation_symbols = []
            KTSSValidator._add_symbol(symbol, self.parser.prefix_map, possible_symbols)
            KTSSValidator._add_symbol(symbol, self.parser.suffix_map, possible_symbols)
            KTSSValidator._add_symbol(
                symbol, self.parser.mutations_map, mutation_symbols
            )

            if KTSSValidator._is_nested(self.parser.mutations_map):
                for key in self.parser.mutations_map:
                    KTSSValidator._add_symbol(
                        symbol, self.parser.mutations_map[key], mutation_symbols
                    )

            possible_symbols = self._filter_possibles(
                current_state, possible_symbols + mutation_symbols
            )

            if len(possible_symbols) < 1:
                if current_state in self.dfa.final_states:
                    return separator.join(result)
                possible_symbols = list(self.dfa.transitions[current_state].keys())

            """TODO: Cambiar cuando se añada el modelo estocástico y usar Viterbi"""
            symbol = possible_symbols[random.randint(0, len(possible_symbols) - 1)]

            result.append(symbol)
            current_state = self.dfa.next_state(symbol, current_state)

        return separator.join(result)

    def _set_mappings(self, parser: ParserVcf):
        """Set the mappings for the prefix, infix, and suffix between an original symbol
        and parsed symbol.

        Parameters
        ----------
        parser: ParserVcf
            Parser where the mappings will be obtained.

        Raise
        -----
        NotImplementedError: When a mapping attribute _inverse_mutations_map is not
        defined on the parser class.
        """
        try:
            self._inverse_mutations_map = parser._inverse_mutations_map()
        except AttributeError:
            raise NotImplementedError(
                "For this metohd is necessary to define prefix_map, mutations_map, suffix_map attributes into the parser class"
            )

    def _string_distances(
        self, string1: str, string2: str, method: Callable = levenshtein
    ) -> Union[int, float]:
        """Mehod to comparing distance between two sequences.

        Parameters
        ----------
        string1: str
            Sequence 1.
        string2: str
            Sequence 2.

        Returns
        -------
        The distance between the two sequences.
        """
        return method(string1, string2)

    @staticmethod
    def _add_symbol(symbol: str, mapping: dict, symbols: list) -> None:
        """If a symbol is present on a mapping as key adds the value to a list. If the
        symbol is not present into the mapping, not add nothing.

        Parameters
        ----------
        symbol: str
            Symbol that will be the key of the mapping.
        mapping: dict
            Map where the value will be obtained.
        symbols: list
            List of symbols to append the value.
        """
        if symbol in mapping:
            symbols.append(mapping[symbol])

    @staticmethod
    def _is_nested(mapping: dict) -> bool:
        """Checks if a map is nested.

        Parameters
        ----------
        mapping: dict
            Map that will be checked.

        Returns
        -------
        True if the map is nested, otherwise false.
        """
        return any(isinstance(i, dict) for i in mapping.values())

    def _filter_possibles(self, state: str, symbols: list) -> list:
        """Filter the symbols that can do a transition in a DFA from a given state.

        Parameters
        ----------
        state: str
            Origin state.
        symbols: list
            List of symbols to filter.

        Returns
        -------
        The filtered list.
        """
        return list(
            filter(lambda symbol: self.dfa.has_transition(state, symbol), symbols)
        )
