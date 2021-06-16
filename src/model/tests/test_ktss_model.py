from unittest import TestCase

from src.model.ktssModel import KTSSModel


class TestKTSSModel(TestCase):
    def setUp(self) -> None:
        self.model = KTSSModel()
        return super().setUp()

    def test__state_in_list(self):
        state = (0, 0, 0)
        lst = [(0, 0, 0), (1, 0, 0), (2, 0, 0), (3, 0, 0)]

        result = self.model._state_in_list(state, lst)

        self.assertTrue(result)

    def test_state_not_in_list(self):
        state = (0, 0, 0)
        lst = [(1, 0, 0), (2, 0, 0), (3, 0, 0)]

        result = self.model._state_in_list(state, lst)

        self.assertFalse(result)

    def test__generate_sigma_k_1(self):
        alphabet = {"C", "G", "T", "A"}
        k = 1
        sigma = {"C", "G", "T", "A"}

        result = self.model._generate_sigma(alphabet, k)

        self.assertEqual(result, sigma)

    def test__generate_sigma_k_2(self):
        alphabet = {"C", "G", "T", "A"}
        k = 2
        sigma = {
            "C",
            "G",
            "T",
            "A",
            "CC",
            "CG",
            "CT",
            "CA",
            "GC",
            "GG",
            "GT",
            "GA",
            "TC",
            "TG",
            "TT",
            "TA",
            "AC",
            "AG",
            "AT",
            "AA",
        }

        result = self.model._generate_sigma(alphabet, k)

        self.assertEqual(result, sigma)

    def test__get_prefix_k_1(self):
        string = "ABCDEFG"
        k = 1
        prefix = ""

        result = self.model._get_prefix(string, k)

        self.assertEqual(result, prefix)

    def test__get_prefix_k_3(self):
        string = "ABCDEFG"
        k = 3
        prefix = "AB"

        result = self.model._get_prefix(string, k)

        self.assertEqual(result, prefix)

    def test__get_prefix_k_greater_string(self):
        string = "ABCDEFG"
        k = 10
        prefix = "ABCDEFG"

        result = self.model._get_prefix(string, k)

        self.assertEqual(result, prefix)

    def test__get_suffix_k_1(self):
        string = "ABCDEFG"
        k = 1
        suffix = ""

        result = self.model._get_suffix(string, k)

        self.assertEqual(result, suffix)

    def test__get_suffix_k_3(self):
        string = "ABCDEFG"
        k = 3
        suffix = "FG"

        result = self.model._get_suffix(string, k)

        self.assertEqual(result, suffix)

    def test__get_suffix_k_greater_string(self):
        string = "ABCDEFG"
        k = 10
        suffix = "ABCDEFG"

        result = self.model._get_suffix(string, k)

        self.assertEqual(result, suffix)

    def test__get_infixes_k_1(self):
        string = "ABCDEFG"
        k = 1
        infixes = [i for i in "ABCDEFG"]

        result = self.model._get_infixes(string, k)

        self.assertEqual(result, infixes)

    def test__get_infixes_k_3(self):
        string = "ABCDEFG"
        k = 3
        infixes = [
            "ABC",
            "BCD",
            "CDE",
            "DEF",
            "EFG",
        ]

        result = self.model._get_infixes(string, k)

        self.assertEqual(result, infixes)

    def test__get_infixes_k_greater_string(self):
        string = "ABCDEFG"
        k = 10
        infixes = ["ABCDEFG"]

        result = self.model._get_infixes(string, k)

        self.assertEqual(result, infixes)

    def test_training_k_2(self):
        samples = ["abba", "aaabba", "bbaaa", "bba"]
        k = 2
        alphabet = {"a", "b"}
        states = {"", "b", "a"}
        transitions = {
            "": {"a": "a", "b": "b"},
            "a": {"a": "a", "b": "b"},
            "b": {"a": "a", "b": "b"},
        }
        initial_state = ""
        final_states = {"a"}
        not_allowed_segments = {"a", "b"}

        result = self.model.training(samples, k, get_not_allowed_segements=True)
        same_transitions = all(map(lambda x: x in transitions, result["transitions"]))

        self.assertEqual(result["alphabet"], alphabet)
        self.assertEqual(result["states"], states)
        self.assertTrue(same_transitions)
        self.assertEqual(result["initial_state"], initial_state)
        self.assertEqual(result["final_states"], final_states)
        self.assertEqual(result["not_allowed_segments"], not_allowed_segments)

    def test_training_k_3(self):
        samples = ["abba", "aaabba", "bbaaa", "bba"]
        k = 3
        alphabet = {"a", "b"}
        states = {"", "ab", "a", "b", "aa", "bb", "ba"}
        transitions = {
            "": {"a": "a", "b": "b"},
            "a": {"b": "ab", "a": "aa"},
            "b": {"b": "bb"},
            "aa": {"b": "ab", "a": "aa"},
            "bb": {"a": "ba"},
            "ba": {"a": "aa"},
            "ab": {"b": "bb"},
        }
        initial_state = ""
        final_states = {"ba", "bba", "aa"}
        not_allowed_segments = {"ab", "a", "b", "aa", "bbb", "bb", "ba", "aba", "bab"}

        result = self.model.training(samples, k, get_not_allowed_segements=True)
        same_transitions = all(map(lambda x: x in transitions, result["transitions"]))

        self.assertEqual(result["alphabet"], alphabet)
        self.assertEqual(result["states"], states)
        self.assertTrue(same_transitions)
        self.assertEqual(result["initial_state"], initial_state)
        self.assertEqual(result["final_states"], final_states)
        self.assertEqual(result["not_allowed_segments"], not_allowed_segments)

    def test__add_transition_empty_transitions(self):
        transitions = {}
        from_state = "a"
        symbol = "b"
        to_state = "c"

        self.model._add_transition(transitions, from_state, symbol, to_state)

        self.assertEqual(transitions, {"a": {"b": "c"}})

    def test__add_transition_not_empty_transitions(self):
        transitions = {"a": {"d": "e"}}
        from_state = "a"
        symbol = "b"
        to_state = "c"

        self.model._add_transition(transitions, from_state, symbol, to_state)

        self.assertEqual(transitions, {"a": {"b": "c", "d": "e"}})

    def test_filter_samples_has_original_False(self):
        sequences = ["a", "b", "a", "b", "a", "b"]
        sequences_obtained = ["a", "b", "a", "b", "a", "b"]

        result = KTSSModel.filter_samples(sequences)

        self.assertEqual(result, sequences_obtained)

    def test_filter_samples_has_original_get_original_False(self):
        sequences = ["a", "b", "a", "b", "a", "b"]
        sequences_obtained = ["a", "b", "a", "b", "a", "b"]

        result = KTSSModel.filter_samples(sequences)

        self.assertEqual(result, sequences_obtained)

    def test_filter_samples_has_original_True_get_original_False(self):
        sequences = ["a", "b", "a", "b", "a", "b"]
        sequences_obtained = ["b", "b", "b"]

        result = KTSSModel.filter_samples(sequences, has_original=True)

        self.assertEqual(result, sequences_obtained)

    def test_filter_samples_has_original_False_get_original_True(self):
        sequences = ["a", "b", "a", "b", "a", "b"]
        sequences_obtained = ["a", "b", "a", "b", "a", "b"]

        result = KTSSModel.filter_samples(sequences)

        self.assertEqual(result, sequences_obtained)

    def test_filter_samples_has_original_True_get_original_True(self):
        sequences = ["a", "b", "a", "b", "a", "b"]
        sequences_obtained = ["a", "a", "a"]

        result = KTSSModel.filter_samples(
            sequences, has_original=True, get_original=True
        )

        self.assertEqual(result, sequences_obtained)
