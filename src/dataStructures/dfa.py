from typing import Union

from sortedcontainers import SortedDict, SortedSet


class DFA(object):
    """
    Data strcuture that represents a deterministic finite automata
    """

    def __init__(
        self,
        states: Union[set, list, tuple],
        alphabet: Union[set, list, tuple],
        transitions: dict,
        initial_state: str,
        final_states: Union[set, list, tuple],
    ):
        self.states = SortedSet(states)
        self.alphabet = SortedSet(alphabet)
        # The transition symbols aren't parsed to a ordered dict because
        # it's not necessary because there is not many transitions per
        # state
        self.transitions = SortedDict(transitions)
        self.initial_state = initial_state
        self.final_states = SortedSet(final_states)

    def add_transition(self, from_state: str, symbol: str, to_state: str):
        """Append a transition into an ordered dict that represent the transitions

        Parameters
        ----------
        from_state: str
            String that represent the source state
        symbol: str
            String that represents the transitions symbol
        to_state: str
            String tht represents the destination state
        """
        if not self.transitions.get(from_state):
            self.transitions[from_state] = SortedDict({})
        self.transitions[from_state][symbol] = to_state

    def next_state(self, symbol: str, state: str) -> str:
        """Parse a symbol from a state and returns the next state

        Parameters
        ----------
        symbol: str
            symbol to be parsed
        state: str
            state where the transition starts

        Raise
        -----
        ValueError: When not exists a transition for given state and symbol

        Returns
        -------
        The next state
        """
        if not self.has_transition(state, symbol):
            raise ValueError(f"The transition for {state} - {symbol} not exists")
        return self.transitions[state][symbol]

    def parse_string(self, string: str, state: bool = False) -> list:
        """Parse a string and returns the sequence of states

        Parameters
        ----------
        string: str
            string to be parsed
        state: bool
            Initial state to generate the sequences

        Raise
        -----
        ValueError: When the string connot be parsed

        Returns
        -------
        Sequence of the followed states
        """
        current_state = state or self.initial_state
        states = [current_state]

        for i in string:
            next_state = self.next_state(i, current_state)
            states.append(next_state)
            current_state = next_state

        return states

    def has_transition(self, state: str, symbol: str) -> bool:
        """Returns true if the transitions from state with symbol exists

        Parameters
        ----------
        state: str
            State of the transition
        symbol: str
            symbol of the transition

        Returns
        -------
        True if the transitions exists otherwise False
        """
        transition = self.transitions.get(state, False)
        return transition and transition.get(symbol, False)
