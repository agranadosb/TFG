import logging
from queue import Queue
from typing import Callable, Union

from sortedcontainers import SortedDict, SortedSet
from src.dataStructures.dfa import DFA
from src.logging.tqdmLoggingHandler import TqdmLoggingHandler
from src.parser.extendedParser import ExtendedParserVcf
from src.parser.parserVcf import ParserVcf
from textdistance import levenshtein
from tqdm import tqdm


class KTSSValidator(object):
    """Operates from a DFA generated from a KTSS model and allows to generate distances
    between original sequences and sequences derived from the original sequences

    Parameters
    ----------
    model: dict
        Dictionary that contains de DFA model
    parser: ParserVCF = ExtendedParserVcf
        Parser of the model
    """

    def __init__(
        self, model: Union[SortedDict, dict], parser: ParserVcf = ExtendedParserVcf
    ):
        self.model = model
        self.infix_symbols = parser.mutations_symbols

        self.set_mappings(parser)

        self.dfa = DFA(
            model["states"],
            model["alphabet"],
            model["transitions"],
            model["initial_state"],
            model["final_states"],
        )

    def set_mappings(self, parser: ParserVcf):
        """Set the mappings for the prefix, infix, and suffix between an original symbol
        and parsed symbol

        Parameters
        ----------
        parser: ParserVcf
            Parser where the mappings will be obtained
        """
        try:
            self.prefix_map = parser.prefix_map
            self.mutations_map = parser.mutations_map
            self.suffix_map = parser.suffix_map
        except AttributeError:
            raise NotImplementedError(
                "For this metohd is necessary to define prefix_map, mutations_map, suffix_map attributes into the parser class"
            )

    def get_next_state(self, sequence: str, state: bool = False) -> str:
        """Parse a sequence and gets the next state after parsing the sequence

        Parameters
        ----------
        sequence: str
            sequence to be parsed

        Raise
        -----
        ValueError: When the sequence connot be parsed

        Returns
        -------
        Next state after parsing the sequence
        """
        states = self.dfa.parse_string(sequence, state=state)

        if states and len(states) > 0:
            return states[len(states) - 1]
        return self.dfa.initial_state

    def generate_sequence(
        self,
        symbol_list: Union[list, tuple],
        state: bool = False,
        loop_tolerance: int = 0,
    ) -> SortedSet:
        """Returns all the possible sequences generated from a given state (or the
        initial if any state is given) and a symbols list.

        For example, if we have the transitions

            - {"a": {"a": "aa", "b": "ab"}}
            - {"aa": {"a": "aaa", "b": "aab"}}

        and the method gets the symbols list ["a"] and the state "a", the method will
        return the sequences:

            - ["a", "aa"]

        Parameters
        ----------
        symbol_list: list, tuple
            List of symbols to generate the sequences
        state: bool = False
            Initial state to generate the sequences
        loop_tolerance: int = 0
            The times that a loop can be executed for each sequence

        Returns
        -------
        List of generated sequences
        """
        state = state or self.dfa.initial_state

        results = []
        steps_remaining = Queue()
        steps_remaining.put(
            {"state": state, "sequence": "", "loop_tolerance": loop_tolerance}
        )

        states_searched = set()
        while not steps_remaining.empty():
            current_step = steps_remaining.get()

            state = current_step["state"]
            sequence = current_step["sequence"]
            tolerance = current_step["loop_tolerance"]

            if state in states_searched:
                if tolerance < 1:
                    continue
                else:
                    tolerance -= 1

            states_searched.add(state)

            for symbol in symbol_list:
                has_transition = self.dfa.has_transition(state, symbol)
                append_sequence = not has_transition and sequence

                if has_transition:
                    next_state = self.dfa.next_state(symbol, state)
                    steps_remaining.put(
                        {
                            "state": next_state,
                            "sequence": sequence + symbol,
                            "loop_tolerance": tolerance,
                        }
                    )

                    transitions = self.dfa.transitions[state].keys()
                    append_sequence = len(transitions) >= 1 and sequence

                if append_sequence:
                    results.append(sequence)

        return SortedSet(results)

    def generate_infixes(
        self, prefix: str, suffix: str, separator: str = "-"
    ) -> SortedSet:
        """Get the sequences generated by the model with a prefix and a suffix and the
        list of infix symbols

        Parameters
        ----------
        prefix: str
            Prefix of the sequence
        suffix str:
            Suffix of the sequence
        separator: str
            Separator between prefix, sequence generated and suffix

        Raise
        -----
        ValueError: When the prefix connot be parsed

        Reutrns
        -------
        List of possible sequences with a given prefix and suffix
        """
        prefix_state = self.get_next_state(prefix)
        sequences = self.generate_sequence(self.infix_symbols, state=prefix_state)

        results = []
        for sequence in sequences:
            state = self.get_next_state(sequence, state=prefix_state)

            if not state:
                continue

            try:
                self.get_next_state(suffix, state=state)
                results.append(f"{prefix}{separator}{sequence}{separator}{suffix}")
            except:
                pass

        return SortedSet(results)

    def string_distances(
        self, string1: str, string2: str, method: Callable = levenshtein
    ) -> Union[int, float]:
        """Mehod to comparing distance between two sequences

        Parameters
        ----------
        string1: str
            Sequence 1
        string2: str
            Sequence 2

        Returns
        -------
        The distance between the two sequences
        """
        return method(string1, string2)

    @staticmethod
    def transform_sequence(sequence: str, mapping: dict) -> str:
        """Transforms a sequence using a dict that repsents a mapping

        Parameters
        ----------
        sequence: str
            Sequence to be mapped
        mapping: dict
            Dictionary that repsents the mapping

        Returns
        -------
        Sequence mapped
        """
        return "".join(map(lambda char: mapping[char], sequence))

    def generate_distances(
        self,
        sequences: Union[list, tuple],
        separator: str = "-",
        prefix_length: int = 5,
        suffix_length: int = 5,
    ) -> SortedDict:
        """Generates all the distances of an infix of a given list of sequences of all
        the possible infixes

        Paameters
        ---------
        sequences: list
            List of sequences
        separator: str
            Separator of the sequence generator

        Returns
        -------
        Dictionary that contains the sequence and its distances
        """
        logger = logging.getLogger()
        tqdm_out = TqdmLoggingHandler(logger, level=logging.INFO)

        result = SortedDict()
        logger.info("Generating validation data")
        for sequence in tqdm(sequences, file=tqdm_out):
            prefix = KTSSValidator.transform_sequence(
                sequence[:prefix_length], self.prefix_map
            )
            suffix = KTSSValidator.transform_sequence(
                sequence[len(sequence) - suffix_length :], self.suffix_map
            )
            infix = KTSSValidator.transform_sequence(
                sequence[prefix_length : len(sequence) - suffix_length],
                self.mutations_map,
            )

            sequence_key = f"{prefix}{separator}{infix}{separator}{suffix}"
            result[sequence_key] = SortedDict()

            try:
                infix_sequences = self.generate_infixes(
                    prefix, suffix, separator=separator
                )
            except:
                result[sequence_key] = False
                continue

            for infix_sequence in infix_sequences:
                infix_string = infix_sequence.split(separator)[1]

                distance = self.string_distances(infix, infix_string)

                result[sequence_key][infix_string] = distance

        return result
