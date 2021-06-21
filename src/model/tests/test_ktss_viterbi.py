from unittest import TestCase

from src.model.ktssViterbi import KTSSViterbi
from src.model.tests.factories import ParserFactoryKTSSValidatorDistances


class TestKTSSValidatorViterbi(TestCase):
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
            "10": {"z": "17"},
            "15": {"z": "16"},
            "16": {"z": "18"},
        }
        probabilities = {
            "1": {"a": 1 / 2, "b": 1 / 2},
            "2": {"a": 1 / 4, "b": 1 / 4, "d": 1 / 4, "l": 1 / 4},
            "3": {"a": 1 / 3, "b": 1 / 3, "d": 1 / 3},
            "4": {"z": 1 / 2, "x": 1 / 2},
            "6": {"z": 1},
            "9": {"z": 1 / 2, "x": 1 / 2},
            "10": {"z": 1},
            "15": {"z": 1},
            "16": {"z": 1},
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
        self.ktss_validator = KTSSViterbi(model)
        self.ktss_validator.parser = ParserFactoryKTSSValidatorDistances
        return super().setUp()

    def test__get_symbol_by_destination(self):
        transitions = {"a": "1", "b": "2"}
        destination = "2"
        origin = "b"

        result = self.ktss_validator._get_symbol_by_destination(
            transitions, destination
        )

        self.assertEqual(result, origin)

    def test__backtrack(self):
        state = "18"
        backward = {0: {"2": "1"}, 1: {"15": "2"}, 2: {"16": "15"}, 3: {"18": "16"}}
        sequence = "alzz"

        result = self.ktss_validator._backtrack(state, backward, 3)

        self.assertEqual(result, sequence)

    def test_annotate_sequence_viterbi(self):
        sequence = "AAAA"
        annotation = "alzz"

        result = self.ktss_validator.annotate_sequence(sequence)

        self.assertEqual(result, annotation)

    def test_annotate_sequence_viterbi_greater(self):
        sequence = "AAAAA"
        annotation = "alzz"

        result = self.ktss_validator.annotate_sequence(sequence)

        self.assertEqual(result, annotation)
