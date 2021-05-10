from unittest import TestCase

from src.model.ktssModel import KTSSModel


class TestKTSSModel(TestCase):
    def setUp(self) -> None:
        self.model = KTSSModel()
        return super().setUp()

    def test_state_in_list(self):
        state = (0, 0, 0)
        lst = [(0, 0, 0), (1, 0, 0), (2, 0, 0), (3, 0, 0)]

        result = self.model.state_in_list(state, lst)

        self.assertTrue(result)

    def test_state_not_in_list(self):
        state = (0, 0, 0)
        lst = [(1, 0, 0), (2, 0, 0), (3, 0, 0)]

        result = self.model.state_in_list(state, lst)

        self.assertFalse(result)

    def test_generate_sigma_k_1(self):
        alphabet = ["C", "G", "T", "A"]
        k = 1
        sigma = ["C", "G", "T", "A"]

        result = self.model.generate_sigma(alphabet, k)

        self.assertEqual(result, sigma)

    def test_generate_sigma_k_2(self):
        alphabet = ["C", "G", "T", "A"]
        k = 2
        sigma = [
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
        ]

        result = self.model.generate_sigma(alphabet, k)

        self.assertEqual(result, sigma)

    def test_get_prefix_k_1(self):
        string = "ABCDEFG"
        k = 1
        prefix = ""

        result = self.model.get_prefix(string, k)

        self.assertEqual(result, prefix)

    def test_get_prefix_k_3(self):
        string = "ABCDEFG"
        k = 3
        prefix = "AB"

        result = self.model.get_prefix(string, k)

        self.assertEqual(result, prefix)

    def test_get_prefix_k_greater_string(self):
        string = "ABCDEFG"
        k = 10
        prefix = "ABCDEFG"

        result = self.model.get_prefix(string, k)

        self.assertEqual(result, prefix)

    def test_get_suffix_k_1(self):
        string = "ABCDEFG"
        k = 1
        suffix = ""

        result = self.model.get_suffix(string, k)

        self.assertEqual(result, suffix)

    def test_get_suffix_k_3(self):
        string = "ABCDEFG"
        k = 3
        suffix = "FG"

        result = self.model.get_suffix(string, k)

        self.assertEqual(result, suffix)

    def test_get_suffix_k_greater_string(self):
        string = "ABCDEFG"
        k = 10
        suffix = "ABCDEFG"

        result = self.model.get_suffix(string, k)

        self.assertEqual(result, suffix)

    def test_get_infixes_k_1(self):
        string = "ABCDEFG"
        k = 1
        infixes = [i for i in "ABCDEFG"]

        result = self.model.get_infixes(string, k)

        self.assertEqual(result, infixes)

    def test_get_infixes_k_3(self):
        string = "ABCDEFG"
        k = 3
        infixes = [
            "ABC",
            "BCD",
            "CDE",
            "DEF",
            "EFG",
        ]

        result = self.model.get_infixes(string, k)

        self.assertEqual(result, infixes)

    def test_get_infixes_k_greater_string(self):
        string = "ABCDEFG"
        k = 10
        infixes = ["ABCDEFG"]

        result = self.model.get_infixes(string, k)

        self.assertEqual(result, infixes)

    def test_training_k_2(self):
        samples = ["abba", "aaabba", "bbaaa", "bba"]
        k = 2
        alphabet = {"a", "b"}
        states = {"", "b", "a"}
        transitions = [
            ["", "a", "a"],
            ["", "b", "b"],
            ["a", "a", "a"],
            ["b", "a", "a"],
            ["b", "b", "b"],
            ["a", "b", "b"],
        ]
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
        transitions = [
            ["", "a", "a"],
            ["a", "b", "ab"],
            ["", "b", "b"],
            ["b", "b", "bb"],
            ["a", "a", "aa"],
            ["aa", "b", "ab"],
            ["bb", "a", "ba"],
            ["ba", "a", "aa"],
            ["ab", "b", "bb"],
            ["aa", "a", "aa"],
        ]
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
