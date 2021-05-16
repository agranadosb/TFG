from unittest import TestCase

from sortedcontainers import SortedSet
from src.model.ktssValidation import KtssValidator


class TestKtssValidator(TestCase):
    def setUp(self) -> None:
        alphabet = {"a", "b"}
        states = {"", "a", "b", "aa", "bb", "ab", "ba"}
        transitions = {
            "": {"a": "a", "b": "b"},
            "a": {"a": "aa", "b": "ab"},
            "b": {"a": "ba", "b": "bb"},
            "aa": {"a": "aaa", "b": "aab"},
            "ab": {"a": "aba", "b": "abb"},
            "ba": {"a": "baa", "b": "bab"},
            "bb": {"a": "bba", "b": "bbb"},
        }
        initial_state = ""
        final_states = {"bb", "aa"}

        model = {
            "alphabet": alphabet,
            "states": states,
            "transitions": transitions,
            "initial_state": initial_state,
            "final_states": final_states,
        }

        self.ktss_validator = KtssValidator(model)

        return super().setUp()

    def test_get_next_state_from_sequence(self):
        prefix = "aa"
        state = "aa"

        result = self.ktss_validator.get_next_state_from_sequence(prefix)

        self.assertEqual(result, state)

    def test_get_next_state_from_empty_string(self):
        prefix = ""
        state = ""

        result = self.ktss_validator.get_next_state_from_sequence(prefix)

        self.assertEqual(result, state)

    def test_generate_sequence_from_list(self):
        symbol_list = ["a", "b"]
        sequences = SortedSet(
            [
                "a",
                "b",
                "aa",
                "ab",
                "ba",
                "bb",
                "aaa",
                "aab",
                "aba",
                "abb",
                "baa",
                "bab",
                "bba",
                "bbb",
            ]
        )

        result = self.ktss_validator.generate_sequence_from_list(symbol_list)

        self.assertEqual(result, sequences)

    def test_generate_sequence_from_list_from_state(self):
        symbol_list = ["a", "b"]
        state = "a"
        sequences = SortedSet(["a", "b", "aa", "ab", "ba", "bb"])

        result = self.ktss_validator.generate_sequence_from_list(
            symbol_list, state=state
        )

        self.assertEqual(result, sequences)

    def test_generate_sequence_from_list_without_symbols(self):
        symbol_list = []
        state = "a"
        sequences = SortedSet()

        result = self.ktss_validator.generate_sequence_from_list(
            symbol_list, state=state
        )

        self.assertEqual(result, sequences)

    def test_generate_sequence_from_list_with_invalid_state(self):
        symbol_list = ["a", "b"]
        state = "c"
        sequences = SortedSet()

        result = self.ktss_validator.generate_sequence_from_list(
            symbol_list, state=state
        )

        self.assertEqual(result, sequences)

    def test_generate_sequence_from_list_with_invalid_symbols(self):
        symbol_list = ["c"]
        sequences = SortedSet()

        result = self.ktss_validator.generate_sequence_from_list(symbol_list)

        self.assertEqual(result, sequences)

    def test_generate_sequence_from_list_with_valid_symbols(self):
        symbol_list = ["a"]
        sequences = SortedSet(["a", "aa", "aaa"])

        result = self.ktss_validator.generate_sequence_from_list(symbol_list)

        self.assertEqual(result, sequences)

    def test_generate_sequence_from_list_example_case(self):
        self.ktss_validator.dfa.transitions = {
            "a": {"a": "aa", "b": "ab"},
            "aa": {"a": "aaa", "b": "aab"},
        }
        symbol_list = ["a"]
        state = "a"
        sequences = SortedSet(["a", "aa"])

        result = self.ktss_validator.generate_sequence_from_list(
            symbol_list, state=state
        )

        self.assertEqual(result, sequences)

    def test_get_model_sequence(self):
        prefix = "a"
        suffix = "b"
        self.ktss_validator.infix_symbols = ["a"]
        sequences = SortedSet(["aab"])

        result = self.ktss_validator.get_model_generated_sequence(prefix, suffix)

        self.assertEqual(result, sequences)

    def test_get_model_sequence_empty_prefix(self):
        prefix = ""
        suffix = "b"
        self.ktss_validator.infix_symbols = ["a"]
        sequences = SortedSet(["ab", "aab"])

        result = self.ktss_validator.get_model_generated_sequence(prefix, suffix)

        self.assertEqual(result, sequences)

    def test_get_model_sequence_empty_prefix_bb_suffix(self):
        prefix = ""
        suffix = "bb"
        self.ktss_validator.infix_symbols = ["a"]
        sequences = SortedSet(["abb"])

        result = self.ktss_validator.get_model_generated_sequence(prefix, suffix)

        self.assertEqual(result, sequences)

    def test_get_model_sequence_empty_suffix(self):
        prefix = "a"
        suffix = ""
        self.ktss_validator.infix_symbols = ["a"]
        sequences = SortedSet(["aa", "aaa"])

        result = self.ktss_validator.get_model_generated_sequence(prefix, suffix)

        self.assertEqual(result, sequences)

    def test_get_model_sequence_all_big_sequences(self):
        self.ktss_validator.dfa.alphabet = {"a", "b", "c"}
        self.ktss_validator.dfa.add_transition("", "c", "c")
        self.ktss_validator.dfa.add_transition("c", "a", "ca")
        self.ktss_validator.dfa.add_transition("c", "b", "cb")
        self.ktss_validator.dfa.add_transition("ca", "a", "caa")
        self.ktss_validator.dfa.add_transition("cb", "a", "cba")
        self.ktss_validator.dfa.add_transition("ca", "b", "cab")
        self.ktss_validator.dfa.add_transition("cb", "b", "cbb")
        self.ktss_validator.dfa.add_transition("caa", "a", "caaa")
        self.ktss_validator.dfa.add_transition("cba", "a", "cbaa")
        self.ktss_validator.dfa.add_transition("cab", "a", "caba")
        self.ktss_validator.dfa.add_transition("cbb", "a", "cbba")
        self.ktss_validator.dfa.add_transition("caa", "b", "caab")
        self.ktss_validator.dfa.add_transition("cba", "b", "cbab")
        self.ktss_validator.dfa.add_transition("cab", "b", "cabb")
        self.ktss_validator.dfa.add_transition("cbb", "b", "cbbb")
        self.ktss_validator.dfa.add_transition("caaa", "c", "caaac")
        self.ktss_validator.dfa.add_transition("cbaa", "c", "cbaac")
        self.ktss_validator.dfa.add_transition("caba", "c", "cabac")
        self.ktss_validator.dfa.add_transition("cbba", "c", "cbbac")
        self.ktss_validator.dfa.add_transition("caab", "c", "caabc")
        self.ktss_validator.dfa.add_transition("cbab", "c", "cbabc")
        self.ktss_validator.dfa.add_transition("cabb", "c", "cabbc")
        self.ktss_validator.dfa.add_transition("cbbb", "c", "cbbbc")
        prefix = "c"
        suffix = "c"
        self.ktss_validator.infix_symbols = ["a", "b"]
        sequences = SortedSet(
            ["caaac", "cbaac", "cabac", "cbbac", "caabc", "cbabc", "cabbc", "cbbbc"]
        )

        result = self.ktss_validator.get_model_generated_sequence(prefix, suffix)

        self.assertEqual(result, sequences)
