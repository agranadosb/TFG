# -*- coding: utf-8 -*-

from unittest import TestCase

from src.utils.arithmetic import _lcm, lcm


class TestArithmetic(TestCase):
    def test__lcm(self):
        number_1 = 3
        number_2 = 10
        lcm = 30

        result = _lcm(number_1, number_2)

        self.assertEqual(result, lcm)

    def test_lcm(self):
        numbers = [3, 5, 7, 11]
        lcm_result = 3 * 5 * 7 * 11

        result = lcm(numbers)

        self.assertEqual(result, lcm_result)
