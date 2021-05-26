from src.parser.extendedParser import ExtendedParserVcf


class ParserFactory(ExtendedParserVcf):
    prefix_map = {"A": "a", "C": "c", "B": "b"}
    mutations_map = {"A": "a", "C": "c", "B": "b"}
    suffix_map = {"A": "a", "C": "c", "B": "b"}

    @classmethod
    def inverse_mutations_map(cls):
        return {}

    def __init__(self):
        pass


class InvalidParserFactory(object):
    def __init__(self):
        pass
