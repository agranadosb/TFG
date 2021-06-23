# -*- coding: utf-8 -*-

import pathlib
from unittest import TestCase

from src.fasta.chromosome import Chromosome


class TestChromosomeCHR1(TestCase):
    def setUp(self) -> None:
        self.static_dir = f"{str(pathlib.Path(__file__).parent.absolute())}/static/"
        self.sequence = "GCATGCATGCATGCATGCATGCATGCATG"
        self.chromosome = Chromosome(
            open(f"{self.static_dir}test.fa", "r"),
            name="chr1",
            line_length=12,
            label_length=6,
            index_start=0,
            length=29,
        )

        return super().setUp()

    def test___getitem__(self):
        sequence = self.sequence[0]

        data = self.chromosome[0]

        self.assertEqual(data, sequence)

    def test___getitem___no_start(self):
        sequence = self.sequence[:24]

        data = self.chromosome[:24]

        self.assertEqual(data, sequence)

    def test___getitem___no_end(self):
        sequence = self.sequence[4:]

        data = self.chromosome[4:]

        self.assertEqual(data, sequence)

    def test___getitem___all(self):
        sequence = self.sequence[:]

        data = self.chromosome[:]

        self.assertEqual(data, sequence)

    def test__get_prefix_at_start(self):
        pos = 5
        length = 7
        prefix = self.sequence[0:5]

        data = self.chromosome[pos:-length]

        self.assertEqual(data, prefix)

    def test__get_prefix(self):
        pos = 13
        length = 7
        prefix = self.sequence[pos - length : pos]

        data = self.chromosome[pos:-length]

        self.assertEqual(data, prefix)

    def test__get_prefix_at_end(self):
        pos = 28
        length = 7
        prefix = self.sequence[pos - length : pos]

        data = self.chromosome[pos:-length]

        self.assertEqual(data, prefix)

    def test__get_nucleotide_index_pos_greater_length(self):
        pos = 1000
        error = False

        try:
            self.chromosome._get_nucleotide_index(pos)
        except IndexError:
            error = True

        self.assertTrue(error)

    def test__get_nucleotide_index_pos_lower_0(self):
        pos = -1
        error = False

        try:
            self.chromosome._get_nucleotide_index(pos)
        except IndexError:
            error = True

        self.assertTrue(error)

    def test__get_nucleotide_index(self):
        pos = 24
        index = 32

        result = self.chromosome._get_nucleotide_index(pos)

        self.assertEqual(result, index)

    def test__get_from_interval(self):
        index = 6
        length = 29
        sequence = "GCATGCATGCATGCATGCATGCATGCATG"

        result = self.chromosome._get_from_interval(index, length)

        self.assertEqual(result, sequence)


class TestChromosomeCHR2(TestCase):
    def setUp(self) -> None:
        self.static_dir = f"{str(pathlib.Path(__file__).parent.absolute())}/static/"
        self.sequence = "TGACTGACTGACTGACTGACTGACTGACTGACTGAC"
        self.chromosome = Chromosome(
            open(f"{self.static_dir}test.fa", "r"),
            name="chr2",
            line_length=12,
            label_length=6,
            index_start=38,
            length=36,
        )

        return super().setUp()

    def test___getitem__(self):
        sequence = self.sequence[0]

        data = self.chromosome[0]

        self.assertEqual(data, sequence)

    def test___getitem___no_start(self):
        sequence = self.sequence[:24]

        data = self.chromosome[:24]

        self.assertEqual(data, sequence)

    def test___getitem___no_end(self):
        sequence = self.sequence[4:]

        data = self.chromosome[4:]

        self.assertEqual(data, sequence)

    def test___getitem___all(self):
        sequence = self.sequence[:]

        data = self.chromosome[:]

        self.assertEqual(data, sequence)

    def test__get_prefix_at_start(self):
        pos = 5
        length = 7
        prefix = self.sequence[0:5]

        data = self.chromosome[pos:-length]

        self.assertEqual(data, prefix)

    def test__get_prefix(self):
        pos = 13
        length = 7
        prefix = self.sequence[pos - length : pos]

        data = self.chromosome[pos:-length]

        self.assertEqual(data, prefix)

    def test__get_prefix_at_end(self):
        pos = 28
        length = 7
        prefix = self.sequence[pos - length : pos]

        data = self.chromosome[pos:-length]

        self.assertEqual(data, prefix)


class TestChromosomeCHR3(TestCase):
    def setUp(self) -> None:
        self.static_dir = f"{str(pathlib.Path(__file__).parent.absolute())}/static/"
        self.sequence = "AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"
        self.chromosome = Chromosome(
            open(f"{self.static_dir}test.fa", "r"),
            name="chr3",
            line_length=12,
            label_length=6,
            index_start=83,
            length=33,
        )

        return super().setUp()

    def test___getitem__(self):
        sequence = self.sequence[0]

        data = self.chromosome[0]

        self.assertEqual(data, sequence)

    def test___getitem___no_start(self):
        sequence = self.sequence[:24]

        data = self.chromosome[:24]

        self.assertEqual(data, sequence)

    def test___getitem___no_end(self):
        sequence = self.sequence[4:]

        data = self.chromosome[4:]

        self.assertEqual(data, sequence)

    def test___getitem___all(self):
        sequence = self.sequence[:]

        data = self.chromosome[:]

        self.assertEqual(data, sequence)

    def test__get_prefix_at_start(self):
        pos = 5
        length = 7
        prefix = self.sequence[0:5]

        data = self.chromosome[pos:-length]

        self.assertEqual(data, prefix)

    def test__get_prefix(self):
        pos = 13
        length = 7
        prefix = self.sequence[pos - length : pos]

        data = self.chromosome[pos:-length]

        self.assertEqual(data, prefix)

    def test__get_prefix_at_end(self):
        pos = 28
        length = 7
        prefix = self.sequence[pos - length : pos]

        data = self.chromosome[pos:-length]

        self.assertEqual(data, prefix)

    def test_get_sequence(self):
        nucleotide = "GC"
        pos = 9
        from_nuc = 9
        to_nuc = 33
        sequence = ["AGCTAGCTA", "GC", "TAGCTAGCTAGCTAGCTAGCTA"]

        result = self.chromosome.sequence(nucleotide, pos, from_nuc, to_nuc)

        self.assertEqual(result, sequence)
