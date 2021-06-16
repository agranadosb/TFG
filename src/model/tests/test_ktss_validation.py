from unittest import TestCase

from sortedcontainers import SortedDict, SortedSet
from src.model.ktssValidation import KTSSValidator
from src.model.tests.factories import InvalidParserFactory, ParserFactory


class TestKTSSValidator(TestCase):
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

        self.ktss_validator = KTSSValidator(model)

        return super().setUp()

    def test__get_next_state(self):
        prefix = "aa"
        state = "aa"

        result = self.ktss_validator._get_next_state(prefix)

        self.assertEqual(result, state)

    def test__get_next_state_from_empty_string(self):
        prefix = ""
        state = ""

        result = self.ktss_validator._get_next_state(prefix)

        self.assertEqual(result, state)

    def test_generate_sequence(self):
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

        result = self.ktss_validator.generate_sequence(symbol_list)

        self.assertEqual(result, sequences)

    def test_generate_sequence_from_state(self):
        symbol_list = ["a", "b"]
        state = "a"
        sequences = SortedSet(["a", "b", "aa", "ab", "ba", "bb"])

        result = self.ktss_validator.generate_sequence(symbol_list, state=state)

        self.assertEqual(result, sequences)

    def test_generate_sequence_without_symbols(self):
        symbol_list = []
        state = "a"
        sequences = SortedSet()

        result = self.ktss_validator.generate_sequence(symbol_list, state=state)

        self.assertEqual(result, sequences)

    def test_generate_sequence_with_invalid_state(self):
        symbol_list = ["a", "b"]
        state = "c"
        sequences = SortedSet()

        result = self.ktss_validator.generate_sequence(symbol_list, state=state)

        self.assertEqual(result, sequences)

    def test_generate_sequence_with_invalid_symbols(self):
        symbol_list = ["c"]
        sequences = SortedSet()

        result = self.ktss_validator.generate_sequence(symbol_list)

        self.assertEqual(result, sequences)

    def test_generate_sequence_with_valid_symbols(self):
        symbol_list = ["a"]
        sequences = SortedSet(["a", "aa", "aaa"])

        result = self.ktss_validator.generate_sequence(symbol_list)

        self.assertEqual(result, sequences)

    def test_generate_sequence_example_case(self):
        self.ktss_validator.dfa.transitions = {
            "a": {"a": "aa", "b": "ab"},
            "aa": {"a": "aaa", "b": "aab"},
        }
        symbol_list = ["a"]
        state = "a"
        sequences = SortedSet(["a", "aa"])

        result = self.ktss_validator.generate_sequence(symbol_list, state=state)

        self.assertEqual(result, sequences)

    def test_get_model_sequence(self):
        prefix = "a"
        suffix = "b"
        self.ktss_validator.infix_symbols = ["a"]
        sequences = SortedSet(["a-a-b"])

        result = self.ktss_validator._generate_infixes(prefix, suffix)

        self.assertEqual(result, sequences)

    def test_get_model_sequence_empty_prefix(self):
        prefix = ""
        suffix = "b"
        self.ktss_validator.infix_symbols = ["a"]
        sequences = SortedSet(["-a-b", "-aa-b"])

        result = self.ktss_validator._generate_infixes(prefix, suffix)

        self.assertEqual(result, sequences)

    def test_get_model_sequence_empty_prefix_bb_suffix(self):
        prefix = ""
        suffix = "bb"
        self.ktss_validator.infix_symbols = ["a"]
        sequences = SortedSet(["-a-bb"])

        result = self.ktss_validator._generate_infixes(prefix, suffix)

        self.assertEqual(result, sequences)

    def test_get_model_sequence_empty_suffix(self):
        prefix = "a"
        suffix = ""
        self.ktss_validator.infix_symbols = ["a"]
        sequences = SortedSet(["a-a-", "a-aa-"])

        result = self.ktss_validator._generate_infixes(prefix, suffix)

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
            [
                "c-aaa-c",
                "c-baa-c",
                "c-aba-c",
                "c-bba-c",
                "c-aab-c",
                "c-bab-c",
                "c-abb-c",
                "c-bbb-c",
            ]
        )

        result = self.ktss_validator._generate_infixes(prefix, suffix)

        self.assertEqual(result, sequences)

    def test_generate_distances_invalid(self):
        self.ktss_validator.infix_symbols = ["a", "b"]
        self.ktss_validator.prefix_map = {"c": "c"}
        self.ktss_validator.mutations_map = {"a": "a", "c": "c"}
        self.ktss_validator.suffix_map = {"c": "c"}
        sequences = [("caaac|a"), ("ccacc|a")]
        distances = SortedDict(
            {
                "caaac": SortedDict({"false-aaa": 0.5}),
                "ccacc": SortedDict({"false-sas": 0.5}),
            }
        )
        prefix_length = 1
        suffix_length = 1

        result = self.ktss_validator.generate_distances(
            sequences, prefix_length=prefix_length, suffix_length=suffix_length
        )

        self.assertEqual(result, distances)

    def test_generate_distances(self):
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
        self.ktss_validator.infix_symbols = ["a", "b"]
        self.ktss_validator.prefix_map = {"a": "a", "c": "c", "b": "b"}
        self.ktss_validator.mutations_map = {"a": "a", "c": "c", "b": "b"}
        self.ktss_validator._inverse_mutations_map = {"a": "a", "c": "c", "b": "b"}
        self.ktss_validator.suffix_map = {"a": "a", "c": "c", "b": "b"}
        self.ktss_validator.parser = ParserFactory
        sequences = [
            "caaac|a",
            "ccacc|a",
            "cabac|a",
            "cbbac|a",
            "ccaac|a",
            "cbabc|a",
            "ccccc|a",
            "caccc|a",
        ]
        prefix_length = 1
        suffix_length = 1
        distances = SortedDict(
            {
                "caaac": SortedDict(
                    {
                        "aaa": 0.0,
                        "aab": 0.007936507936507936,
                        "aba": 0.007936507936507936,
                        "abb": 0.015873015873015872,
                        "baa": 0.007936507936507936,
                        "bab": 0.015873015873015872,
                        "bba": 0.015873015873015872,
                        "bbb": 0.023809523809523808,
                    }
                ),
                "cabac": SortedDict(
                    {
                        "aaa": 0.007936507936507936,
                        "aab": 0.015873015873015872,
                        "aba": 0.0,
                        "abb": 0.007936507936507936,
                        "baa": 0.015873015873015872,
                        "bab": 0.015873015873015872,
                        "bba": 0.007936507936507936,
                        "bbb": 0.015873015873015872,
                    }
                ),
                "caccc": SortedDict(
                    {
                        "aaa": 0.015873015873015872,
                        "aab": 0.015873015873015872,
                        "aba": 0.015873015873015872,
                        "abb": 0.015873015873015872,
                        "baa": 0.023809523809523808,
                        "bab": 0.023809523809523808,
                        "bba": 0.023809523809523808,
                        "bbb": 0.023809523809523808,
                    }
                ),
                "cbabc": SortedDict(
                    {
                        "aaa": 0.015873015873015872,
                        "aab": 0.007936507936507936,
                        "aba": 0.015873015873015872,
                        "abb": 0.015873015873015872,
                        "baa": 0.007936507936507936,
                        "bab": 0.0,
                        "bba": 0.015873015873015872,
                        "bbb": 0.007936507936507936,
                    }
                ),
                "cbbac": SortedDict(
                    {
                        "aaa": 0.015873015873015872,
                        "aab": 0.023809523809523808,
                        "aba": 0.007936507936507936,
                        "abb": 0.015873015873015872,
                        "baa": 0.007936507936507936,
                        "bab": 0.015873015873015872,
                        "bba": 0.0,
                        "bbb": 0.007936507936507936,
                    }
                ),
                "ccaac": SortedDict(
                    {
                        "aaa": 0.007936507936507936,
                        "aab": 0.015873015873015872,
                        "aba": 0.015873015873015872,
                        "abb": 0.023809523809523808,
                        "baa": 0.007936507936507936,
                        "bab": 0.015873015873015872,
                        "bba": 0.015873015873015872,
                        "bbb": 0.023809523809523808,
                    }
                ),
                "ccacc": SortedDict(
                    {
                        "aaa": 0.015873015873015872,
                        "aab": 0.015873015873015872,
                        "aba": 0.023809523809523808,
                        "abb": 0.023809523809523808,
                        "baa": 0.015873015873015872,
                        "bab": 0.015873015873015872,
                        "bba": 0.023809523809523808,
                        "bbb": 0.023809523809523808,
                    }
                ),
                "ccccc": SortedDict(
                    {
                        "aaa": 0.023809523809523808,
                        "aab": 0.023809523809523808,
                        "aba": 0.023809523809523808,
                        "abb": 0.023809523809523808,
                        "baa": 0.023809523809523808,
                        "bab": 0.023809523809523808,
                        "bba": 0.023809523809523808,
                        "bbb": 0.023809523809523808,
                    }
                ),
            }
        )

        result = self.ktss_validator.generate_distances(
            sequences, prefix_length=prefix_length, suffix_length=suffix_length
        )

        self.assertEqual(result, distances)

    def test__set_mappings_valid(self):
        parser = ParserFactory()

        is_valid = False
        try:
            self.ktss_validator._set_mappings(parser)
            is_valid = True
        except NotImplementedError:
            pass

        self.assertTrue(is_valid)

    def test__set_mappings_invalid(self):
        is_invalid = False
        try:
            self.ktss_validator._set_mappings(InvalidParserFactory)
        except NotImplementedError:
            is_invalid = True

        self.assertTrue(is_invalid)

    def test__get_minimum_distances(self):
        distances = {"a": {"a": 1, "b": 0, "c": 2}, "b": {}, "c": False}
        minimum_distances = {"a": {"b": 0}, "b": -1, "c": -1}

        result = KTSSValidator._get_minimum_distances(distances)

        self.assertEqual(result, minimum_distances)
