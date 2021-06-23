from src.model.ktssValidation import KTSSValidator


class KTSSViterbi(KTSSValidator):
    def annotate_sequence(self, sequence: str, separator: str = "") -> str:
        """Gets a string sequence, and annotates it using Viterbi algorthm.

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

        We can get two results:

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
        forward_porb = {self.dfa.initial_state: 1}
        backward = {}
        leafs = set(self.dfa.initial_state)
        bad_leafs = set()

        count = 0
        for symbol in sequence:
            auxiliar_leafs = set()
            level_backward = {}
            for current_state in leafs:
                possible_symbols = self._get_possible_symbols(current_state, symbol)
                if not possible_symbols:
                    bad_leafs.add(current_state)
                    continue

                for possible_symbol in possible_symbols:
                    next_state = self.dfa.transitions[current_state][possible_symbol]
                    probability = self.dfa.probabilities[current_state][possible_symbol]
                    forward_porb[next_state] = forward_porb[current_state] * probability
                    level_backward[next_state] = current_state

                    auxiliar_leafs.add(next_state)
            if len(level_backward.keys()) > 0:
                backward[count] = level_backward
            leafs = auxiliar_leafs
            count += 1

        if len(leafs) < 1:
            leafs = bad_leafs

        states_probabilities = {state: forward_porb[state] for state in leafs}
        state = max(states_probabilities, key=states_probabilities.get)

        return self._backtrack(
            state, backward, len(backward.keys()) - 1, separator=separator
        )

    def _get_symbol_by_destination(self, transitions: dict, destination: str) -> str:
        """Gets a symbol transition from a transition by de destination state.

        For instnce, if the transition is `{"a": "1", "b": "2"}` and the destination is
        `"2"` the result is `"b"`.

        Parameters
        ----------
        transitions: dict
            Transitions.
        destination: str
            One transition destination.

        Returns
        -------
        The symbol that goes to the destination.
        """
        return list(transitions.keys())[list(transitions.values()).index(destination)]

    def _backtrack(
        self,
        state: str,
        backward: dict,
        length: int,
        separator: str = "",
    ) -> str:
        """Backtracking through a dictionary to get the sequence of generated symbols
        using the states.

        For example, if the backward dictionary is:

        ```python
            {0: {"2": "1"}, 1: {"15": "2"}, 2: {"16": "15"}, 3: {"18": "16"}}
        ```

        The state is `"18"` and the dfa transitions are:

        ```python
            "1": {"a": "2", "b": "3"},
            "2": {"a": "4", "b": "5", "d": "6", "l": "15"},
            "3": {"a": "7", "b": "8", "d": "9"},
            "4": {"z": "10", "x": "11"},
            "6": {"z": "12"},
            "9": {"z": "13", "x": "14"},
            "10": {"z": "17"},
            "15": {"z": "16"},
            "16": {"z": "18"},
        ```

        The method will return `"alzz"`.

        Parameters
        ----------
        state: str
            State where de backtracking will start.
        backward: dict
            Backward dictionary.
        length: int
            Length of the backward dictionary.
        separator: str = ""
            Separator between the symbols of the annotated sequence.

        Returns
        -------
        The generated sequence.
        """
        res = []
        for i in range(length, -1, -1):
            if not state in backward[i]:
                continue
            previous_state = backward[i][state]
            symbol = self._get_symbol_by_destination(
                self.dfa.transitions[previous_state], state
            )
            res = [symbol] + res
            state = previous_state
        return separator.join(res)
