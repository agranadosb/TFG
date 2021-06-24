from src.dataStructures.dfaStochastic import DFAStochastic


class WatsonCrickAutomata(object):
    """Class that represents a [Watson-Crick Automata](https://arxiv.org/pdf/1602.05721.pdf) 1-limited"""

    def __init__(
        self,
        alphabet,
        pairs,
        states,
        initial_state,
        final_states,
        transitions,
        probabilities,
    ) -> None:
        self.alphabet = alphabet
        self.pairs = pairs
        self.states = states
        self.initial_state = initial_state
        self.final_states = final_states
        self.transitions = transitions
        self.probabilities = probabilities
        super().__init__()

    def _hash_pair(self, pair: tuple, symbol: str = "-") -> str:
        """Creates a unic string by a pair.

        Parameters
        ----------
        pair: tuple
            Pair to be hashed.
        symbol: str = "-"
            Symbol to divide the hash.

        Returns
        -------
        Return the hashed pair.
        """
        return f"{pair[0]}{symbol}{pair[1]}"

    def add_transition(
        self, origin: str, pair: tuple, destination: str, probability: bool = False
    ) -> dict:
        """Append a transition into an ordered dict that represent the transitions.

        Parameters
        ----------
        from_state: str
            String that represent the source state.
        pair: str
            Pair symbol to represent the transition.
        to_state: str
            String tht represents the destination state.
        probability: bool = False
            Probability of the transitition.

        Returns
        -------
        Dictionary with the transition.
        """
        self.alphabet.add(pair[0])
        self.alphabet.add(pair[1])
        self.pairs.append(pair)
        self.states.add(origin)
        self.states.add(destination)

        if not self.transitions.get(origin, False):
            self.transitions[origin] = {}
            self.probabilities[origin] = {}

        hashed_pair = self._hash_pair(pair)
        transition = {hashed_pair: destination}

        self.transitions[origin][hashed_pair] = destination
        if not probability:
            DFAStochastic._add_probability(self.probabilities[origin], hashed_pair)
        else:
            self.probabilities[origin][hashed_pair] = probability

        return transition

    def next_state(self, pair: tuple, origin: str) -> str:
        """Parse a symbol from a state and returns the next state.

        Parameters
        ----------
        pair: tuple
            Pair.
        state: str
            State where the transition starts.

        Returns
        -------
        The next state.
        """
        return self.transitions[origin][self._hash_pair(pair)]

    def parse_dfa(self, dfa, mappings):
        mappings_keys = mappings.keys()
        for origin in dfa.transitions:
            transition = dfa.transitions[origin]
            available_symbols = list(
                filter(lambda x: x in mappings_keys, dfa.transitions[origin].keys())
            )
            for symbol in available_symbols:
                self.add_transition(
                    origin,
                    (symbol, mappings[symbol]),
                    transition[symbol],
                    probability=dfa.probabilities[origin][symbol],
                )

        return self.transitions

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

        Parameters
        ----------
        sequence: str
            Sequence to be annotated.
        separator: str = ""
            Separator between the symbols of the annotated sequence.

        Returns
        -------
        Annotated sequence
        """
        result = []
        current_state = self.initial_state
        for symbol in sequence:
            symbols = list(
                map(lambda x: x.split("-"), self.transitions[current_state].keys())
            )

            symbols = [self._hash_pair(i) for i in symbols if i[1].upper() == symbol]
            if not symbols:
                return separator.join(result)

            state_probabilities = self.probabilities[current_state]

            probabilities = [
                (possible_symbol, state_probabilities[possible_symbol])
                for possible_symbol in symbols
            ]
            max_probability_symbol = max(probabilities, key=lambda item: item[1])[0]
            new_symbol = max_probability_symbol.split("-")

            result.append(new_symbol[0])

            current_state = self.next_state(new_symbol, current_state)

        return separator.join(result)
