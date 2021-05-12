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
        
        Returns
        -------
        The next state
        """
        transition = self.transitions.get(state, False)
        if not transition or not transition.get(symbol, False):
            raise ValueError(f"The transition for {state} - {symbol} not exists")
        return transition[symbol]
