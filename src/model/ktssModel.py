# -*- coding: utf-8 -*-

import functools
import itertools
import json
import operator

from src.model.abstractModel import AbstractModel


class KTSSModel(AbstractModel):
    def __init__(self, save_path=False, restore_path=False):
        super().__init__(save_path=save_path, restore_path=restore_path)
        self.model = False
        self.trainer_name = "ktt"

    def state_in_list(self, state, lst):
        """Checks if a state is in a list

        Parameters
        ----------
        state: tuple
            State to be searched (s0, s1, s2)
        lst: list
            list where the state will be searched

        Returns
        -------
        True if the state is in the list, false either
        """
        for state_list in lst:
            if (
                state_list[0] == state[0]
                and state_list[1] == state[1]
                and state_list[2] == state[2]
            ):
                return True
        return False

    def generate_sigma(self, alphabet, k):
        """Generates sigma for a k given, for example, if k = 0 and alphabet = [A, B]
        it generates [A, B, AA, AB, BA, BB]

        Parameters
        ----------
        alphabet: list
            alphabet to generate sigma
        k: int
            length of the maximum strings of sigma

        Returns
        -------
        sigma
        """
        aux = 1
        words = [] + alphabet
        while aux < k:
            words += [pair[0] + pair[1] for pair in itertools.product(words, alphabet)]
            aux += 1

        return words

    def get_prefix(self, string, k):
        """Returns the prefix of a string

        Parameters
        ----------
        string: str
            Strong where the prefix will be getted
        k: int
            Length of the prefix - 1

        Returns
        -------
        The prefix of the string
        """
        return string[: k - 1]

    def get_suffix(self, string, k):
        """Returns the suffix of a string

        Parameters
        ----------
        string: str
            String where the suffix will be getted
        k: int
            Length of the suffix - 1

        Returns
        -------
        The suffix of the string
        """
        if k >= len(string):
            return string
        return string[len(string) - k + 1 :]

    def get_infixes(self, string, k):
        """Generates all the possible infixes of a string for a lengt k given

        Parameters
        ----------
        string: str
            String where the infixes will be getted
        k: int
            Length of the infixes

        Rreturns
        --------
        List of all the infixes
        """
        if k >= len(string):
            return [string]
        ini = 0
        fin = k
        result = []
        while fin <= len(string):
            result.append(string[ini:fin])
            ini += 1
            fin += 1
        return result

    def training(self, samples, k):
        """Generates a ktss model from the samples and a k given

        Parameters
        ----------
        samples: list
            List of samples to training the model
        k: int
            Parameter of the ktss model

        Returns
        -------
        {
            "states": states,
            "alphabet": alphabet,
            "transitions": transitions,
            "initial_state": initial_state,
            "final_states": final_states,
            "not_allowed_segments": not_allowed_segments,
        }
        """
        alphabet = list(set(functools.reduce(operator.add, samples)))

        greater_or_equal_than_k = []
        lower_than_k = []
        for sample in samples:
            if len(sample) >= k:
                greater_or_equal_than_k.append(sample)
            else:
                lower_than_k.append(sample)

        prefixes = [] + lower_than_k
        suffixes = [] + lower_than_k
        infixes = []
        for sample in greater_or_equal_than_k:
            prefixes.append(self.get_prefix(sample, k))
            suffixes.append(self.get_suffix(sample, k))
            infixes += self.get_infixes(sample, k)

        prefixes = set(prefixes)
        suffixes = set(suffixes)
        infixes = set(infixes)

        sigma_k = self.generate_sigma(alphabet, k)
        not_allowed_segments = set(sigma_k) - set(infixes)

        q = [""]
        s = []
        q0 = ""

        for prefix in prefixes:
            s.append(["", prefix[0], prefix[0]])
            for char_index in range(len(prefix)):
                q.append([prefix[: char_index + 1]])
                s.append(
                    [prefix[:char_index], prefix[char_index], prefix[: char_index + 1]]
                )

        for infix in infixes:
            q.append([infix[: k - 1], infix[2:k]])
            s.append([infix[: k - 1], infix[k - 1], infix[1:k]])

        q = set(
            list(
                map(
                    lambda x: x if type(x) == str else x[0],
                    filter(lambda x: type(x) == str or all(x), q),
                )
            )
        )

        s_aux = []
        for i in s:
            if not self.state_in_list(i, s_aux):
                s_aux.append(i)

        s = s_aux

        self.model = {
            "states": q,
            "alphabet": set(alphabet),
            "transitions": s,
            "initial_state": q0,
            "final_states": suffixes,
            "not_allowed_segments": not_allowed_segments,
        }

        return self.model

    def parser(self):
        return self.parser_class

    def trainer(self):
        return self.training

    def model(self):
        return self.model

    def saver(self):
        if not self.save_path:
            raise AttributeError("Saver path is not defined")

        if not self.model:
            raise ValueError("The model is not trained")

        with open(f"{self.save_path}ktss-model.json", "w") as outfile:
            model_for_json = {
                "states": list(self.model["states"]),
                "alphabet": list(self.model["alphabet"]),
                "transitions": self.model["transitions"],
                "initial_state": self.model["initial_state"],
                "final_states": list(self.model["final_states"]),
                "not_allowed_segments": list(self.model["not_allowed_segments"]),
            }
            json.dump(model_for_json, outfile)

    def loader(self):
        if not self.restore_path:
            raise AttributeError("Loader path is not defined")

        with open(self.save_path) as json_file:
            self.model = json.load(json_file)
