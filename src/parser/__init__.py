# -*- coding: utf-8 -*-

"""
**ParserVcf**
-------------

Parses data from the VCF and FASTA files and prepares that data for a machine
learning model.

The parser needs to define a method named `method`, which generates a new sequence
from the original one.

For example, if we have the sequence `abbc` (a, b, c) and the mutation `aa` this
`method` changes it into other representation, for example to `smms` (s, mm, s)
where the mutation will be added (and transformed) in the middle of a tuple, being
the start and the end of that tuple the prefix and the suffix respectively:

```python
    ("s", "mm", "s") = method(("a", "b", "c"), "bb")
```

The parser will have static methods to send and retrieve the parsed and transformed
sequence to a file in a specified format. To do this, the parser will need to be defined the methods:

- ** sequence_to_string **: Transforms a sequence (in a tuple format) to a string
format, adding to this string the original sequence (in a string format), the
mutation and a prefix.
- ** retrive_sequence **: Gets a sequence in a string format and returns it in a
tuple format. This method will get the string sequences that had been transformed by
`sequence_to_string` method.
- ** retrive_string_sequence **: Gets a sequence in a string format and returns it
in a different string format. This method will get the string sequences that had
been transformed by `sequence_to_string` method.

**ExtendedParser** *extends ParserVcf*
--------------------------------------

Parses data from the files VCF and FASTA and prepares that data for a machine
learning model in a file.

The extended parser gets a sequence in a tuple format (prefix, infix, suffix) and
maps each character of the sequence to another depending if the character is in the
prefix, infix or suffix:

- Prefix mapping:
    - A -> q
    - C -> w
    - G -> e
    - T -> r
- Infix mapping:
    - A -> a
    - C -> s
    - G -> d
    - T -> f
- Suffix mapping:
    - A -> z
    - C -> x
    - G -> c
    - T -> v

So, for instance, for a given sequence ("ACGT", "ACGT", "ACGT") the parsers change
it to:

- ("qwer", "asdf", "zxcv")

The transformed sequence is named **extended** sequence.

**MutationParser** *extends ExtendedParser*
-------------------------------------------

Parses data from the files VCF and FASTA and prepares that data for a machine
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
"""
from . import extendedParser, parserVcf

__pdoc__ = {}

__pdoc__["tests"] = False
