import logging
from queue import Queue

from sortedcontainers import SortedDict, SortedSet
from src.dataStructures.dfa import DFA
from src.logging.tqdmLoggingHandler import TqdmLoggingHandler
from src.parser.extendedParser import ExtendedParserVcf
from textdistance import levenshtein
from tqdm import tqdm


class KTSSValidator(object):
    """Operates from a DFA generated from a KTSS model and allows to generate distances
    between original sequences and sequences derived from the original sequences

    Parameters
    ----------
    model: dict
        Dictionary that contains de DFA model
    [parser: ParserVCF = ExtendedParserVcf]
        Parser of the model
    """

    def __init__(self, model, parser=ExtendedParserVcf):
        self.model = model
        self.infix_symbols = parser.mutations_symbols
        self.dfa = DFA(
            model["states"],
            model["alphabet"],
            model["transitions"],
            model["initial_state"],
            model["final_states"],
        )

    def get_next_state(self, sequence, state=False):
        """Parse a sequence and gets the next state after parsing the sequence

        Parameters
        ----------
        string: str
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

    def generate_sequence(self, symbol_list, state=False):
        """Returns all the possible sequences generated from a given state (or the
        initial if any state is given) and a symbols list.

        For example, if we have the transitions
            - {"a": {"a": "aa", "b": "ab"}}
            - {"aa": {"a": "aaa", "b": "aab"}}
        and the method gets the symbols list ["a"] and the state "a", the method will
        return the sequences:
            - ["a", "aa"].

        Parameters
        ----------
        symbol_list: list
            List of symbols to generate the sequences
        [state: bool]
            Initial state to generate the sequences

        Returns
        -------
        List of generated sequences
        """
        state = state or self.dfa.initial_state

        results = []
        steps_remaining = Queue()
        steps_remaining.put({"state": state, "sequence": ""})

        while not steps_remaining.empty():
            current_step = steps_remaining.get()

            state = current_step["state"]
            sequence = current_step["sequence"]

            for symbol in symbol_list:
                has_transition = self.dfa.has_transition(state, symbol)
                append_sequence = not has_transition and sequence

                if has_transition:
                    next_state = self.dfa.next_state(symbol, state)
                    steps_remaining.put(
                        {"state": next_state, "sequence": sequence + symbol}
                    )

                    transitions = self.dfa.transitions[state].keys()
                    append_sequence = len(transitions) >= 1 and sequence

                if append_sequence:
                    results.append(sequence)

        return SortedSet(results)

    def generate_infixes(self, prefix, suffix, separator="-"):
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

    def string_distances(self, string1, string2, method=levenshtein):
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

    def generate_distances(self, sequences, separator="-"):
        """Generates all the distances of an infix of a given list of sequences of all
        the possible infixes

        Paameters
        ---------
        sequences: list
            List of sequences
        [separator: str]
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
            prefix = sequence[0]
            infix = sequence[1]
            suffix = sequence[2]
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
