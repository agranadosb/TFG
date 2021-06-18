from unittest import TestCase

from src.dataStructures.dfaStochastic import DFAStochastic


class TestDFAStochastic(TestCase):
    def setUp(self) -> None:
        alphabet = {"a", "b"}
        states = {"", "a", "b", "aa", "bb", "ab", "ba"}
        transitions = {}
        probabilities = {}
        initial_state = ""
        final_states = {"bb", "aa"}

        self.dfa = DFAStochastic(
            states, alphabet, transitions, initial_state, final_states, probabilities
        )
        return super().setUp()

    def test__add_probability(self):
        transition = {"a": 1 / 4, "b": 1 / 4, "c": 1 / 2}
        symbol = "a"
        new_transition = {"a": 2 / 5, "b": 1 / 5, "c": 2 / 5}

        DFAStochastic._add_probability(transition, symbol)

        self.assertEqual(transition, new_transition)

    def test__add_probability_empty(self):
        transition = {}
        symbol = "a"
        new_transition = {"a": 1}

        DFAStochastic._add_probability(transition, symbol)

        self.assertEqual(transition, new_transition)

    def test__add_probability_not_present(self):
        transition = {"b": 1 / 2, "c": 1 / 2}
        symbol = "a"
        new_transition = {"a": 1 / 3, "b": 1 / 3, "c": 1 / 3}

        DFAStochastic._add_probability(transition, symbol)

        self.assertEqual(transition, new_transition)

    def test_add_transition_empty_transitions(self):
        from_state = "a"
        symbol = "b"
        to_state = "ab"
        transition = {"a": {"b": 1}}

        self.dfa.add_transition(from_state, symbol, to_state)

        self.assertEqual(self.dfa.probabilities, transition)

    def test_add_transition_not_empty_transitions(self):
        self.dfa.add_transition("a", "a", "aa")
        from_state = "a"
        symbol = "b"
        to_state = "ab"
        transition = {"a": {"a": 1 / 2, "b": 1 / 2}}

        self.dfa.add_transition(from_state, symbol, to_state)

        self.assertEqual(self.dfa.probabilities, transition)
