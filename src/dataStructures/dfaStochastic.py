from fractions import Fraction
from typing import OrderedDict, Union

from sortedcontainers import SortedDict
from src.dataStructures.dfa import DFA
from src.utils.arithmetic import lcm


class DFAStochastic(DFA):
    def __init__(
        self,
        states: Union[set, list, tuple],
        alphabet: Union[set, list, tuple],
        transitions: dict,
        initial_state: str,
        final_states: Union[set, list, tuple],
        probabilities: dict,
    ):
        super().__init__(states, alphabet, transitions, initial_state, final_states)

        self.probabilities = SortedDict(probabilities)

    @staticmethod
    def _add_probability(
        transition: Union[OrderedDict, dict], symbol: str
    ) -> OrderedDict:
        """Update the probabilities of a transition when a new symbol is added.

        For instance, if our transition is:

        ```python
            {"a": 1 / 4, "b": 1 / 4, "c": 1 / 2}
        ```

        And our symbol is `a` the result is:

        ```python
            {"a": 2 / 5, "b": 1 / 5, "c": 2 / 5}
        ```

        Parameters
        ----------
        transition: OrderedDict, dict
            Transition to be updated.
        symbol: str
            Symbol that was added.
        """
        if not transition.get(symbol, False) and not len(transition.values()):
            transition[symbol] = 1
            return transition

        if not transition.get(symbol, False):
            transition[symbol] = 0

        lcm_value = lcm(
            [Fraction(number).denominator for number in transition.values()]
        )

        for transition_symbol in transition:
            fraction = Fraction(transition[transition_symbol])

            add = 0
            if transition_symbol == symbol:
                add = 1

            new_numerator = (
                lcm_value / fraction.denominator
            ) * fraction.numerator + add
            new_denominator = lcm_value + 1

            transition[transition_symbol] = new_numerator / new_denominator

        return transition

    def add_transition(self, from_state: str, symbol: str, to_state: str):
        """Append a transition into an ordered dict that represent the transitions.

        Parameters
        ----------
        from_state: str
            String that represent the source state.
        symbol: str
            String that represents the transitions symbol.
        to_state: str
            String tht represents the destination state.
        """
        if not self.transitions.get(from_state):
            self.transitions[from_state] = SortedDict({})
            self.probabilities[from_state] = SortedDict({})

        self.transitions[from_state][symbol] = to_state
        DFAStochastic._add_probability(self.probabilities[from_state], symbol)
