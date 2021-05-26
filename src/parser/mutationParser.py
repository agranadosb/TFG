from typing import Union

import Levenshtein
from src.parser.extendedParser import ExtendedParserVcf
from src.utils.genomics import generate_dict_values

MUTATIONS_INSERTION_SYMBOLS = ["t", "y", "u", "i"]
MUTATIONS_EARSED_SYMBOLS = ["g", "h", "j", "k"]
MUTATIONS_SUBSITUTION_SYMBOLS = ["a", "s", "d", "f"]

MUTATION_SYMBOLS = (
    MUTATIONS_INSERTION_SYMBOLS
    + MUTATIONS_EARSED_SYMBOLS
    + MUTATIONS_SUBSITUTION_SYMBOLS
)

INSERT = "insert"
DELETE = "delete"
REPLACE = "replace"

MUTATION_TYPES = [
    (INSERT, MUTATIONS_INSERTION_SYMBOLS),
    (DELETE, MUTATIONS_EARSED_SYMBOLS),
    (REPLACE, MUTATIONS_SUBSITUTION_SYMBOLS),
]


class MutationParser(ExtendedParserVcf):
    """Parser that transforms sequences into mutation type sequences and operates
    between them.

    A mutation type sequence is a sequence that has a different symbol for each type of
    mutation.

    It means that we can get three types of mutations:

        Insertion
        Erase
        Substitution

    And, if we have A,C,G and T symbols, each symbol has three possible mutations:

        A -> [A_insertion, A_earsed, A_substitution]
        C -> [C_insertion, C_earsed, C_substitution]
        G -> [G_insertion, G_earsed, G_substitution]
        T -> [T_insertion, T_earsed, T_substitution]

    In this case, this method differentiates between prefix symbols and suffix symbols,
    so the transformation of each symbol in each position is:

    - **Prefix**
        - A -> q
        - C -> w
        - G -> e
        - T -> r
    - **Infix**
        - **Insertion**
            - A -> t
            - C -> y
            - G -> u
            - T -> i
        - **Erase**
            - A -> g
            - C -> h
            - G -> j
            - T -> k
        - **Substitution**
            - A -> a
            - C -> s
            - G -> d
            - T -> f
    - **Suffix**
        - A -> z
        - C -> x
        - G -> c
        - T -> v
    """

    name: str = "mutation-type"

    mutations_symbols: Union[list, tuple] = (
        MUTATIONS_INSERTION_SYMBOLS
        + MUTATIONS_EARSED_SYMBOLS
        + MUTATIONS_SUBSITUTION_SYMBOLS
    )

    mutations_map: dict = {
        operation: generate_dict_values(symbols)
        for operation, symbols in MUTATION_TYPES
    }

    @classmethod
    def inverse_mutations_map(cls):
        res = {}
        for operation in cls.mutations_map:
            for symbol in cls.mutations_map[operation]:
                res[cls.mutations_map[operation][symbol]] = symbol
        return res

    @classmethod
    def method(cls, sequence: Union[tuple, list], mutation: str) -> tuple:
        """Generates the mutation type sequence from a sequence:

        Parameters
        ----------
        sequence : tuple
            Sequence to simplify
        mutation : str
            Mutation sequence

        Returns
        -------
        Sequence in a mutation type way
        """
        infix = sequence[1]

        operations = Levenshtein.editops(infix, mutation)

        mutation_type_sequence = []
        for operation in operations:
            mutation_type = operation[0]

            reference_sequence = infix
            index = 1
            if mutation_type == INSERT:
                reference_sequence = mutation
                index = 2

            reference_symbol = reference_sequence[operation[index]]

            mutation_type_symbol = cls.mutations_map[mutation_type][reference_symbol]

            mutation_type_sequence.append(mutation_type_symbol)

        left = [cls.prefix_map[nucletid.upper()] for nucletid in sequence[0]]
        right = [cls.suffix_map[nucletid.upper()] for nucletid in sequence[2]]

        return (left, mutation_type_sequence, right)
