from sortedcontainers import SortedDict, SortedSet


class DFA(object):
    """
    Data strcuture dat represents a deterministic finite automata
    """
    def __init__(
        self,
        states,
        alphabet,
        transitions,
        initial_state,
        final_states
    ):
        self.states = SortedSet(states)
        self.alphabet = SortedSet(alphabet)
        self.transitions = SortedDict(transitions)
        self.initial_state = initial_state
        self.final_states = SortedSet(final_states)

    def append_to_transitions(self, from_state, symbol, to_state):
        """Append a transition into an ordered dict that represent the transitions

        Parameters
        ----------
        from_state: str
            String that represent the source state
        symbol: str
            String that represents the transitions symbol
        to_state
            String tht represents the destination state
        """
        if not self.transitions.get(from_state):
            self.transitions[from_state] = SortedDict({})
        self.transitions[from_state][symbol] = to_state

    def next_state(self, symbol, state):
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
        transition = self.transitions.get(state, False)
        if not transition or not transition.get(symbol, False):
            raise ValueError(f"The transition for {state} - {symbol} not exists")
        return transition[symbol]

    def parse_string(self, string):
        """Parse a string and returns the sequence of states

        Parameters
        ----------
        string: str
            string to be parsed
        
        Raise
        -----
        ValueError: When the string connot be parsed

        Returns
        -------
        Sequence of the followed states
        """
        current_state = self.initial_state
        states = [current_state]

        for i in string:
            next_state = self.next_state(i, current_state)
            states.append(next_state)
            current_state = next_state
        
        return states
