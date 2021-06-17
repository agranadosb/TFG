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
    """Parses data from the files VCF and FASTA and prepares that data for a machine
    learning model in a file.

    The **mutation parser** is similar to the **extended parser**, the difference is the
    symbols present in the infix (or mutation) of the transformed sequence. In the
    mutation parser case, each type of modification between the reference and the
    mutation has a different symbol.

    This means that we can get three types of mutations:

    - **Insertion**
    - **Erase**
    - **Substitution**

    So, if we have "A", "C", "G" and "T" symbols for the infix, each symbol has three
    possible mutations:

        A -> ["A_insertion", "A_earsed", "A_substitution"]
        C -> ["C_insertion", "C_earsed", "C_substitution"]
        G -> ["G_insertion", "G_earsed", "G_substitution"]
        T -> ["T_insertion", "T_earsed", "T_substitution"]

    In this case, this method differentiates between prefix symbols and suffix symbols
    (as the extended parser), so the transformation of each symbol in each position is:

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

    The transformed sequence is named the **mutation** sequence.

    Parameters
    ----------
    vcf_path: str
        Path of the vcf file.
    fasta_path: str
        Path of the fasta file.
    """

    name: str = "mutation-type"

    mutations_symbols: Union[list, tuple] = (
        MUTATIONS_INSERTION_SYMBOLS
        + MUTATIONS_EARSED_SYMBOLS
        + MUTATIONS_SUBSITUTION_SYMBOLS
    )
    """List of mutation symbols."""

    mutations_map: dict = {
        operation: generate_dict_values(symbols)
        for operation, symbols in MUTATION_TYPES
    }
    """ Mapping between nucleotides and the mutation symbols divided in each operation:

        ```python
            {
                operation_1: {...} # mapping
                operation_2: {...} # mapping
                ...
            }
        ```
    """

    @classmethod
    def _inverse_mutations_map(cls):
        """Inverse of `mutations_map`"""
        res = {}
        for operation in cls.mutations_map:
            for symbol in cls.mutations_map[operation]:
                res[cls.mutations_map[operation][symbol]] = symbol
        return res

    @classmethod
    def method(cls, sequence: Union[tuple, list], mutation: str) -> tuple:
        """Parse the sequence with a given mutation to a new sequence with different
        noation, in this case a **mutation** sequence.

        For instance, if the sequence is

        ```python
            ("ACGT","ACGT","ACGT")
        ```

        and our mutation is:

        ```python
            "GTTCAC"
        ```

        the method changes it to:

        ```python
            (
                ["q", "w", "e", "r"],
                ["u", "a", "d", "t", "f"],
                ["z", "x", "c", "v"]
            )
        ```

        Parameters
        ----------
        sequence : tuple
            Sequence.
        mutation : str
            Mutation sequence.

        Returns
        -------
        Transformed sequence.
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
