from unittest import TestCase, load_tests

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

    def test_append_to_transitions_empty_transitions(self):
        from_state = "a"
        symbol = "b"
        to_state = "ab"

        self.dfa.append_to_transitions(from_state, symbol, to_state)

        self.assertEqual(self.dfa.transitions, {"a": {"b": "ab"}})

    def test_append_to_transitions_not_empty_transitions(self):
        self.dfa.append_to_transitions("a", "a", "aa")
        from_state = "a"
        symbol = "b"
        to_state = "ab"

        self.dfa.append_to_transitions(from_state, symbol, to_state)

        self.assertEqual(self.dfa.transitions, {"a": {"a": "aa", "b": "ab"}})

    def test_next_state_valid(self):
        self.dfa.append_to_transitions("a", "a", "aa")
        from_state = "a"
        symbol = "a"
        next_state = "aa"

        result = self.dfa.next_state(symbol, from_state)

        self.assertEqual(result, next_state)

    def test_next_state_invalid(self):
        self.dfa.append_to_transitions("a", "a", "aa")
        from_state = "a"
        symbol = "b"

        invalid = False
        try:
            self.dfa.next_state(symbol, from_state)
        except ValueError:
            invalid = True

        self.assertTrue(invalid)

    def test_parse_string_valid(self):
        self.dfa.append_to_transitions("", "a", "a")
        self.dfa.append_to_transitions("a", "b", "ab")
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
