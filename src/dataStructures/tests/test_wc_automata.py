from unittest import TestCase

from src.dataStructures.dfaStochastic import DFAStochastic
from src.dataStructures.watsonCrickAutomata import WatsonCrickAutomata


class TestWCA(TestCase):
    def setUp(self) -> None:
        alphabet = set()
        pairs = []
        states = {"1"}
        initial_state = "1"
        final_states = set()
        transitions = {}
        probabilities = {}

        self.wca = WatsonCrickAutomata(
            alphabet,
            pairs,
            states,
            initial_state,
            final_states,
            transitions,
            probabilities,
        )
        return super().setUp()

    def test__hash_pair(self):
        pair = ("a", "b")
        hashed = "a-b"

        result = self.wca._hash_pair(pair)

        self.assertEqual(result, hashed)

    def test_add_transition(self):
        pair = ("a", "b")
        origin = "1"
        destination = "2"
        transition = {f"{pair[0]}-{pair[1]}": destination}

        result = self.wca.add_transition(origin, pair, destination)

        self.assertEqual(result, transition)

    def test_add_transition_multiple(self):
        pair = ("a", "b")
        origin = "1"
        destination = "2"
        origin2 = "2"
        destination2 = "3"
        transition = {f"{pair[0]}-{pair[1]}": destination2}

        result = self.wca.add_transition(origin, pair, destination)
        result = self.wca.add_transition(origin2, pair, destination2)

        self.assertEqual(result, transition)

    def test_next_state(self):
        pair = ("a", "b")
        origin = "1"
        destination = "2"
        self.wca.add_transition(origin, pair, destination)

        result = self.wca.next_state(pair, origin)

        self.assertEqual(result, destination)

    def test_parse_dfa(self):
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
        dfa = DFAStochastic(
            states, alphabet, transitions, initial_state, final_states, probabilities
        )
        mappings = {
            "a": "A",
            "c": "C",
            "b": "B",
            "l": "A",
            "d": "C",
            "z": "A",
            "x": "C",
            "c": "B",
        }
        transitions = {
            "1": {"a-A": "2", "b-B": "3"},
            "15": {"z-A": "16"},
            "2": {"a-A": "4", "b-B": "5", "d-C": "6", "l-A": "15"},
            "3": {"a-A": "7", "b-B": "8", "d-C": "9"},
            "4": {"z-A": "10", "x-C": "11"},
            "6": {"z-A": "12"},
            "9": {"z-A": "13", "x-C": "14"},
        }
        result = self.wca.parse_dfa(dfa, mappings)

        self.assertEqual(result, transitions)


class TestWCAParse(TestCase):
    def setUp(self) -> None:
        alphabet = set()
        pairs = []
        states = {"1"}
        initial_state = "1"
        final_states = set()
        transitions = {}
        probabilities = {}

        self.wca = WatsonCrickAutomata(
            alphabet,
            pairs,
            states,
            initial_state,
            final_states,
            transitions,
            probabilities,
        )

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
        dfa = DFAStochastic(
            states, alphabet, transitions, initial_state, final_states, probabilities
        )
        mappings = {
            "a": "A",
            "c": "C",
            "b": "B",
            "l": "A",
            "d": "C",
            "z": "A",
            "x": "C",
            "c": "B",
        }
        self.wca.parse_dfa(dfa, mappings)
        return super().setUp()

    def test_annotate_sequence_empty(self):
        sequence = ""
        annotated = ""

        result = self.wca.annotate_sequence(sequence)

        self.assertEqual(result, annotated)

    def test_2(self):
        sequence = "A"
        annotated = "a"

        result = self.wca.annotate_sequence(sequence)

        self.assertEqual(result, annotated)

    def test_3(self):
        sequence = "AAA"
        annotated1 = "aaz"
        annotated2 = "alz"

        results = [self.wca.annotate_sequence(sequence) for _ in range(1000)]

        self.assertIn(annotated1, results)
        self.assertNotIn(annotated2, results)

    def test_4(self):
        sequence = "AC"
        annotated = "ad"

        result = self.wca.annotate_sequence(sequence)

        self.assertEqual(result, annotated)

    def test_5(self):
        sequence = "AA"
        annotated_mutation = "al"
        annotated_not_muutation = "aa"

        results = [self.wca.annotate_sequence(sequence) for _ in range(1000)]

        self.assertNotIn(annotated_mutation, results)
        self.assertIn(annotated_not_muutation, results)

    def test_6(self):
        sequence = "AAA"
        annotated_mutation = "alz"
        annotated_not_muutation = "aaz"

        results = [self.wca.annotate_sequence(sequence) for _ in range(1000)]

        self.assertNotIn(annotated_mutation, results)
        self.assertIn(annotated_not_muutation, results)
