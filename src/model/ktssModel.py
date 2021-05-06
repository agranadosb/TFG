# -*- coding: utf-8 -*-

import functools
import itertools
import operator

from src.model.abstractModel import AbstractModel
from src.model.parserVcf import ParserVcf


class KTSSModel(AbstractModel):
    def __init__(self):
        self.model = False
        self.parser = ParserVcf
        self.trainer = 'k-tt'

    def state_in_list(self, state, lst):
        for state_list in lst:
            if (
                state_list[0] == state[0]
                and state_list[1] == state[1]
                and state_list[2] == state[2]
            ):
                return True
        return False

    def generate_sigma(self, alphabet, k):
        aux = 1
        words = [] + alphabet
        while aux < k:
            words += [pair[0] + pair[1] for pair in itertools.product(words, alphabet)]
            aux += 1

        return words

    def get_prefix(self, string, k):
        return string[: k - 1]

    def get_sufix(self, string, k):
        return string[len(string) - k + 1 :]

    def get_infixes(self, string, k):
        ini = 0
        fin = k
        result = []
        while fin <= len(string):
            result.append(string[ini:fin])
            ini += 1
            fin += 1
        return result

    def training(self, samples, k):
        alphabet = list(set(functools.reduce(operator.add, samples)))

        greater_or_equal_than_k = []
        lower_than_k = []
        for sample in samples:
            if len(sample) >= k:
                greater_or_equal_than_k.append(sample)
            else:
                lower_than_k.append(sample)

        prefixes = [] + lower_than_k
        sufixes = [] + lower_than_k
        infixes = []
        for sample in greater_or_equal_than_k:
            prefixes.append(self.get_prefix(sample))
            sufixes.append(self.get_sufix(sample))
            infixes += self.get_infixes(sample)

        prefixes = set(prefixes)
        sufixes = set(sufixes)
        infixes = set(infixes)

        sigma_k = self.generate_sigma()
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
            "final_states": sufixes,
            "not_allowed_segments": not_allowed_segments,
        }
