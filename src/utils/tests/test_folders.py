# -*- coding: utf-8 -*-

from unittest import TestCase

from src.utils.folders import parse_route


class TestFastaReader(TestCase):
    def test_parse_route_with_char(self):
        route = "route"
        parsed_route = "route/"

        result = parse_route(route)

        self.assertEqual(parsed_route, result)
    
    def test_parse_route_false(self):
        route = ""
        parsed_route = ""

        result = parse_route(route)

        self.assertEqual(parsed_route, result)
    
    def test_parse_route_correct(self):
        route = "route/"

        result = parse_route(route)

        self.assertEqual(route, result)
