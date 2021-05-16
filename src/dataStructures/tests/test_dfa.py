from unittest import TestCase

from src.dataStructures.dfa import DFA


class TestDFA(TestCase):
    def setUp(self) -> None:
        alphabet = {"a", "b"}
        states = {"", "a", "b", "aa", "bb", "ab", "ba"}
        transitions = {}
        initial_state = ""
        final_states = {"bb", "aa"}

        self.dfa = DFA(states, alphabet, transitions, initial_state, final_states)
        return super().setUp()

    def test_add_transition_empty_transitions(self):
        from_state = "a"
        symbol = "b"
        to_state = "ab"

        self.dfa.add_transition(from_state, symbol, to_state)

        self.assertEqual(self.dfa.transitions, {"a": {"b": "ab"}})

    def test_add_transition_not_empty_transitions(self):
        self.dfa.add_transition("a", "a", "aa")
        from_state = "a"
        symbol = "b"
        to_state = "ab"

        self.dfa.add_transition(from_state, symbol, to_state)

        self.assertEqual(self.dfa.transitions, {"a": {"a": "aa", "b": "ab"}})

    def test_next_state_valid(self):
        self.dfa.add_transition("a", "a", "aa")
        from_state = "a"
        symbol = "a"
        next_state = "aa"

        result = self.dfa.next_state(symbol, from_state)

        self.assertEqual(result, next_state)

    def test_next_state_invalid(self):
        self.dfa.add_transition("a", "a", "aa")
        from_state = "a"
        symbol = "b"

        invalid = False
        try:
            self.dfa.next_state(symbol, from_state)
        except ValueError:
            invalid = True

        self.assertTrue(invalid)

    def test_parse_string_valid(self):
        self.dfa.add_transition("", "a", "a")
        self.dfa.add_transition("a", "b", "ab")
        string = "ab"

        result = self.dfa.parse_string(string)

        self.assertEqual(result, ["", "a", "ab"])

    def test_parse_string_invalid(self):
        string = "ab"

        invalid = False
        try:
            self.dfa.parse_string(string)
        except ValueError:
            invalid = True

        self.assertTrue(invalid)

    def test_has_transition_valid(self):
        self.dfa.add_transition("", "a", "a")
        state = ""
        symbol = "a"

        result = self.dfa.has_transition(state, symbol)

        self.assertTrue(result)


    def test_has_transition_invalid_state(self):
        self.dfa.add_transition("", "a", "a")
        state = "b"
        symbol = "a"

        result = self.dfa.has_transition(state, symbol)

        self.assertFalse(result)

    def test_has_transition_invalid_symbol(self):
        self.dfa.add_transition("", "a", "a")
        state = ""
        symbol = "b"

        result = self.dfa.has_transition(state, symbol)

        self.assertFalse(result)
