from unittest import TestCase

from src.model.ktssValidation import KTSSValidator
from src.model.tests.factories import (
    InvalidParserFactory, ParserFactory, ParserFactoryKTSSValidatorAnnotate,
    ParserFactoryKTSSValidatorAnnotateNested,
    ParserFactoryKTSSValidatorDistances)


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
        probabilities = {
            "": {"a": 1 / 2, "b": 1 / 2},
            "a": {"a": 1 / 2, "b": 1 / 2},
            "b": {"a": 1 / 2, "b": 1 / 2},
            "aa": {"a": 1 / 2, "b": 1 / 2},
            "ab": {"a": 1 / 2, "b": 1 / 2},
            "ba": {"a": 1 / 2, "b": 1 / 2},
            "bb": {"a": 1 / 2, "b": 1 / 2},
        }
        initial_state = ""
        final_states = {"bb", "aa"}

        model = {
            "alphabet": alphabet,
            "states": states,
            "transitions": transitions,
            "initial_state": initial_state,
            "final_states": final_states,
            "probabilities": probabilities,
        }

        self.ktss_validator = KTSSValidator(model)

        return super().setUp()

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

    def test__string_distances_match(self):
        string1 = "AAA"
        string2 = "AAA"
        difference = 0

        result = self.ktss_validator._string_distances(string1, string2)

        self.assertEqual(result, difference)

    def test__string_distances_not_match(self):
        string1 = "AAA"
        string2 = "AAB"
        difference = 1

        result = self.ktss_validator._string_distances(string1, string2)

        self.assertEqual(result, difference)

    def test__add_symbol_true(self):
        symbol = "a"
        mapping = {"a": "s"}
        symbols = []
        result = ["s"]

        KTSSValidator._add_symbol(symbol, mapping, symbols)

        self.assertEqual(symbols, result)

    def test__add_symbol_false(self):
        symbol = "v"
        mapping = {"a": "s"}
        symbols = []
        result = []

        KTSSValidator._add_symbol(symbol, mapping, symbols)

        self.assertEqual(symbols, result)

    def test__is_nested_true(self):
        mapping = {"a": {}}

        result = KTSSValidator._is_nested(mapping)

        self.assertTrue(result)

    def test__is_nested_false(self):
        mapping = {}

        result = KTSSValidator._is_nested(mapping)

        self.assertFalse(result)

    def test__filter_possibles(self):
        state = ""
        symbols = ["a", "b", "c", "d"]
        possible_symbols = ["a", "b"]

        result = self.ktss_validator._filter_possibles(state, symbols)

        self.assertEqual(result, possible_symbols)

    def test__filter_possibles_empty(self):
        state = ""
        symbols = ["c", "d"]
        possible_symbols = []

        result = self.ktss_validator._filter_possibles(state, symbols)

        self.assertEqual(result, possible_symbols)


class TestKTSSValidatorAnnotate(TestCase):
    def setUp(self) -> None:
        alphabet = {"a", "b", "d"}
        states = {"1", "2", "3", "4", "5", "6", "7"}
        transitions = {
            "1": {"a": "2", "b": "3"},
            "2": {"a": "4", "b": "5", "d": "6", "l": "15"},
            "3": {"a": "7", "b": "8", "d": "9"},
            "4": {"z": "10", "x": "11"},
            "6": {"z": "12"},
            "9": {"z": "13", "x": "14"},
            "15": {"z": "16"},
        }
        probabilities = {
            "1": {"a": 1 / 2, "b": 1 / 2},
            "2": {"a": 1 / 4, "b": 1 / 4, "d": 1 / 4, "l": 1 / 4},
            "3": {"a": 1 / 3, "b": 1 / 3, "d": 1 / 3},
            "4": {"z": 1 / 2, "x": 1 / 2},
            "6": {"z": 1},
            "9": {"z": 1 / 2, "x": 1 / 2},
            "15": {"z": 1},
        }
        initial_state = "1"
        final_states = {}
        model = {
            "alphabet": alphabet,
            "states": states,
            "transitions": transitions,
            "initial_state": initial_state,
            "final_states": final_states,
            "probabilities": probabilities,
        }
        self.ktss_validator = KTSSValidator(model)
        self.ktss_validator.parser = ParserFactoryKTSSValidatorAnnotate
        return super().setUp()

    def test_annotate_sequence_empty(self):
        sequence = ""
        annotated = ""

        result = self.ktss_validator.annotate_sequence(sequence)

        self.assertEqual(result, annotated)

    def test_2(self):
        sequence = "A"
        annotated = "a"

        result = self.ktss_validator.annotate_sequence(sequence)

        self.assertEqual(result, annotated)

    def test_3(self):
        sequence = "AAA"
        annotated1 = "aaz"
        annotated2 = "alz"

        results = [self.ktss_validator.annotate_sequence(sequence) for _ in range(1000)]

        self.assertIn(annotated1, results)
        self.assertIn(annotated2, results)

    def test_4(self):
        sequence = "AC"
        annotated = "ad"

        result = self.ktss_validator.annotate_sequence(sequence)

        self.assertEqual(result, annotated)

    def test_5(self):
        sequence = "AA"
        annotated_mutation = "al"
        annotated_not_muutation = "aa"

        results = [self.ktss_validator.annotate_sequence(sequence) for _ in range(1000)]

        self.assertIn(annotated_mutation, results)
        self.assertIn(annotated_not_muutation, results)

    def test_6(self):
        sequence = "AAA"
        annotated_mutation = "alz"
        annotated_not_muutation = "aaz"

        results = [self.ktss_validator.annotate_sequence(sequence) for _ in range(1000)]

        self.assertIn(annotated_mutation, results)
        self.assertIn(annotated_not_muutation, results)

    def test_9(self):
        sequence = "Z"
        annotated1 = "a"
        annotated2 = "b"

        results = [self.ktss_validator.annotate_sequence(sequence) for _ in range(1000)]

        self.assertIn(annotated1, results)
        self.assertIn(annotated2, results)


class TestKTSSValidatorAnnotateNested(TestCase):
    def setUp(self) -> None:
        alphabet = {"a", "b", "d"}
        states = {"1", "2", "3", "4", "5", "6", "7"}
        transitions = {
            "1": {"a": "2", "b": "3"},
            "2": {"a": "4", "b": "5", "d": "6", "l": "15"},
            "3": {"a": "7", "b": "8", "d": "9"},
            "4": {"z": "10", "x": "11"},
            "6": {"z": "12"},
            "9": {"z": "13", "x": "14"},
            "15": {"z": "16"},
        }
        probabilities = {
            "1": {"a": 1 / 2, "b": 1 / 2},
            "2": {"a": 1 / 4, "b": 1 / 4, "d": 1 / 4, "l": 1 / 4},
            "3": {"a": 1 / 3, "b": 1 / 3, "d": 1 / 3},
            "4": {"z": 1 / 2, "x": 1 / 2},
            "6": {"z": 1},
            "9": {"z": 1 / 2, "x": 1 / 2},
            "15": {"z": 1},
        }
        initial_state = "1"
        final_states = {}
        model = {
            "alphabet": alphabet,
            "states": states,
            "transitions": transitions,
            "initial_state": initial_state,
            "final_states": final_states,
            "probabilities": probabilities,
        }
        self.ktss_validator = KTSSValidator(model)
        self.ktss_validator.parser = ParserFactoryKTSSValidatorAnnotateNested
        return super().setUp()

    def test_7(self):
        sequence = "AA"
        annotated1 = "aa"
        annotated2 = "al"

        results = [self.ktss_validator.annotate_sequence(sequence) for _ in range(1000)]

        self.assertIn(annotated1, results)
        self.assertIn(annotated2, results)


class TestKTSSValidatorDistances(TestCase):
    def setUp(self) -> None:
        alphabet = {"a", "b", "d", "z", "x"}
        states = {"1", "2", "3", "4", "5", "6", "7"}
        transitions = {
            "1": {"a": "2", "b": "3"},
            "2": {"a": "4", "b": "5", "d": "6", "l": "15"},
            "3": {"a": "7", "b": "8", "d": "9"},
            "4": {"z": "10", "x": "11"},
            "6": {"z": "12"},
            "9": {"z": "13", "x": "14"},
            "15": {"z": "16"},
        }
        probabilities = {
            "1": {"a": 1 / 2, "b": 1 / 2},
            "2": {"a": 1 / 4, "b": 1 / 4, "d": 1 / 4, "l": 1 / 4},
            "3": {"a": 1 / 3, "b": 1 / 3, "d": 1 / 3},
            "4": {"z": 1 / 2, "x": 1 / 2},
            "6": {"z": 1},
            "9": {"z": 1 / 2, "x": 1 / 2},
            "15": {"z": 1},
        }
        initial_state = "1"
        final_states = {}
        model = {
            "alphabet": alphabet,
            "states": states,
            "transitions": transitions,
            "initial_state": initial_state,
            "final_states": final_states,
            "probabilities": probabilities
        }
        self.ktss_validator = KTSSValidator(model)
        self.ktss_validator.parser = ParserFactoryKTSSValidatorDistances
        return super().setUp()

    def test_generate_distances(self):
        sequences = [("AAA", "aaz")]
        annotated1 = {"AAA": 0, "error": 0}
        annotated2 = {"AAA": 1, "error": 1 / 3}

        results = [
            self.ktss_validator.generate_distances(sequences) for _ in range(1000)
        ]

        self.assertIn(annotated1, results)
        self.assertIn(annotated2, results)

    def test_generate_distances_empty(self):
        sequences = []
        empty = {}

        result = self.ktss_validator.generate_distances(sequences)

        self.assertEqual(result, empty)
