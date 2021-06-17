from src.parser.extendedParser import ExtendedParserVcf


class ParserFactory(ExtendedParserVcf):
    prefix_map = {"A": "a", "C": "c", "B": "b"}
    mutations_map = {"A": "a", "C": "c", "B": "b"}
    suffix_map = {"A": "a", "C": "c", "B": "b"}

    @classmethod
    def _inverse_mutations_map(cls):
        return {}

    def __init__(self):
        pass


class ParserFactoryKTSSValidatorAnnotate(ExtendedParserVcf):
    prefix_map = {"A": "a", "C": "c", "B": "b"}
    mutations_map = {"A": "l", "C": "d"}
    suffix_map = {"A": "z", "C": "x", "B": "c"}

    mutations_symbols = ["d"]

    @classmethod
    def _inverse_mutations_map(cls):
        return {}

    def __init__(self):
        pass


class ParserFactoryKTSSValidatorAnnotateNested(ExtendedParserVcf):
    prefix_map = {"A": "a", "C": "c", "B": "b"}
    mutations_map = {"p": {"A": "l", "C": "d"}}
    suffix_map = {"A": "z", "C": "x", "B": "c"}

    mutations_symbols = ["d"]

    @classmethod
    def _inverse_mutations_map(cls):
        return {}

    def __init__(self):
        pass


class ParserFactoryKTSSValidatorDistances(ExtendedParserVcf):
    prefix_map = {"A": "a", "C": "c", "B": "b"}
    mutations_map = {"A": "l", "C": "d"}
    suffix_map = {"A": "z", "C": "x", "B": "c"}

    mutations_symbols = ["d"]

    @classmethod
    def _inverse_mutations_map(cls):
        return {"l": "A", "d": "C"}

    def __init__(self):
        pass


class InvalidParserFactory(object):
    def __init__(self):
        pass
