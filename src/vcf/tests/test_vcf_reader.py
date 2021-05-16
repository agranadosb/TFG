# -*- coding: utf-8 -*-

import pathlib
from unittest import TestCase

from src.vcf.vcfReader import VcfMutationsReader


class TestVcfMutationsReader(TestCase):
    def setUp(self) -> None:
        self.static_dir = f"{str(pathlib.Path(__file__).parent.absolute())}/static/"
        self.reader = VcfMutationsReader(
            f"{self.static_dir}vcfTest.vcf",
            f"{self.static_dir}test.fa.gz",
        )

        return super().setUp()

    def test_set_fasta_line_length(self):
        line_length = 50

        self.assertEqual(self.reader.fasta_file_line_length, line_length)

    def test_set_fasta_file(self):
        path = f"{self.static_dir}test.fa"

        self.assertEqual(self.reader.fasta_file.name, path)

    def test_get_chromosome_index(self):
        chromosme_1 = "chr1"
        chromosme_2 = "chr2"
        chromosme_3 = "chr3"
        result_1 = 0
        result_2 = len(f">{chromosme_1}\n") + 51 * 10
        result_3 = len(f">{chromosme_1}\n") + len(f">{chromosme_2}\n") + 51 * 10 * 2

        index_1 = self.reader.get_chromosome_index(chromosme_1)
        index_2 = self.reader.get_chromosome_index(chromosme_2)
        index_3 = self.reader.get_chromosome_index(chromosme_3)

        self.assertEqual(index_1, result_1)
        self.assertEqual(index_2, result_2)
        self.assertEqual(index_3, result_3)

    def test_set_chromosme_start_data_first(self):
        chromosme_label = "chr1"
        label_length = len(f">{chromosme_label}\n")
        line_start = 0
        index_start = 0

        self.assertTrue(self.reader.fatsa_data.get(chromosme_label, True))
        self.assertEqual(
            self.reader.fatsa_data[chromosme_label]["label_length"],
            label_length,
        )
        self.assertEqual(
            self.reader.fatsa_data[chromosme_label]["line_start"],
            line_start,
        )
        self.assertEqual(
            self.reader.fatsa_data[chromosme_label]["index_start"],
            index_start,
        )

    def test_set_chromosme_start_data_middle(self):
        chromosme_label = "chr2"
        label_length = len(f">{chromosme_label}\n")
        line_start = 11
        index_start = 516

        self.assertTrue(self.reader.fatsa_data.get(chromosme_label, True))
        self.assertEqual(
            self.reader.fatsa_data[chromosme_label]["label_length"],
            label_length,
        )
        self.assertEqual(
            self.reader.fatsa_data[chromosme_label]["line_start"],
            line_start,
        )
        self.assertEqual(
            self.reader.fatsa_data[chromosme_label]["index_start"],
            index_start,
        )

    def test_set_chromosme_start_data_last(self):
        chromosme_label = "chr3"
        label_length = len(f">{chromosme_label}\n")
        line_start = 22
        index_start = 1032

        self.assertTrue(self.reader.fatsa_data.get(chromosme_label, True))
        self.assertEqual(
            self.reader.fatsa_data[chromosme_label]["label_length"],
            label_length,
        )
        self.assertEqual(
            self.reader.fatsa_data[chromosme_label]["line_start"],
            line_start,
        )
        self.assertEqual(
            self.reader.fatsa_data[chromosme_label]["index_start"],
            index_start,
        )

    def test_set_chromosme_end_data_first(self):
        chromosme_label = "chr1"
        index_ends = 516
        next_label_length = len(f">chr2\n")
        file_ends = False
        index = 0
        chromosme_length = 50 * 10

        self.assertTrue(self.reader.fatsa_data.get(chromosme_label, True))
        self.assertEqual(
            self.reader.fatsa_data[chromosme_label]["index_ends"],
            index_ends,
        )
        self.assertEqual(
            self.reader.fatsa_data[chromosme_label]["next_label_length"],
            next_label_length,
        )
        self.assertEqual(
            self.reader.fatsa_data[chromosme_label]["file_ends"],
            file_ends,
        )
        self.assertEqual(
            self.reader.fatsa_data[chromosme_label]["index"],
            index,
        )
        self.assertEqual(
            self.reader.fatsa_data[chromosme_label]["chromosme_length"],
            chromosme_length,
        )

    def test_set_chromosme_end_data_last(self):
        chromosme_label = "chr3"
        index_ends = 1548
        next_label_length = 0
        file_ends = True
        index = 2
        chromosme_length = 50 * 10

        self.assertTrue(self.reader.fatsa_data.get(chromosme_label, True))
        self.assertEqual(
            self.reader.fatsa_data[chromosme_label]["index_ends"],
            index_ends,
        )
        self.assertEqual(
            self.reader.fatsa_data[chromosme_label]["next_label_length"],
            next_label_length,
        )
        self.assertEqual(
            self.reader.fatsa_data[chromosme_label]["file_ends"],
            file_ends,
        )
        self.assertEqual(
            self.reader.fatsa_data[chromosme_label]["index"],
            index,
        )
        self.assertEqual(
            self.reader.fatsa_data[chromosme_label]["chromosme_length"],
            chromosme_length,
        )

    def test_get_nucleotid_index_lower_0(self):
        chromosme = "chr1"
        pos = -1

        invalid = False
        try:
            self.reader.get_nucleotid(pos, chromosme)
        except IndexError:
            invalid = True

        self.assertTrue(invalid)

    def test_get_nucleotid_index_greater_chromosme_length(self):
        chromosme = "chr1"
        pos = 500

        invalid = False
        try:
            self.reader.get_nucleotid(pos, chromosme)
        except IndexError:
            invalid = True

        self.assertTrue(invalid)

    def test_get_nucleotid_first_line_first(self):
        chromosme = "chr1"
        pos = 0
        index = len(f">{chromosme}\n") + pos
        nucleotid = "T"

        index_result = self.reader.get_nucleotid(pos, chromosme)
        self.reader.fasta_file.seek(index_result, 0)
        nucleotid_result = self.reader.fasta_file.read(1)

        self.assertEqual(index, index_result)
        self.assertEqual(nucleotid, nucleotid_result)

    def test_get_nucleotid_first_line_second(self):
        chromosme = "chr1"
        pos = 1
        index = len(f">{chromosme}\n") + pos
        nucleotid = "A"

        index_result = self.reader.get_nucleotid(pos, chromosme)
        self.reader.fasta_file.seek(index_result, 0)
        nucleotid_result = self.reader.fasta_file.read(1)

        self.assertEqual(index, index_result)
        self.assertEqual(nucleotid, nucleotid_result)

    def test_get_nucleotid_first_line_last(self):
        chromosme = "chr1"
        pos = 49
        index = len(f">{chromosme}\n") + pos
        nucleotid = "C"

        index_result = self.reader.get_nucleotid(pos, chromosme)
        self.reader.fasta_file.seek(index_result, 0)
        nucleotid_result = self.reader.fasta_file.read(1)

        self.assertEqual(index, index_result)
        self.assertEqual(nucleotid, nucleotid_result)

    def test_get_nucleotid_second_line_first(self):
        chromosme = "chr1"
        pos = 50
        index = len(f">{chromosme}\n") + pos + 1
        nucleotid = "A"

        index_result = self.reader.get_nucleotid(pos, chromosme)
        self.reader.fasta_file.seek(index_result, 0)
        nucleotid_result = self.reader.fasta_file.read(1)

        self.assertEqual(index, index_result)
        self.assertEqual(nucleotid, nucleotid_result)

    def test_get_nucleotid_second_line_second(self):
        chromosme = "chr1"
        pos = 51
        index = len(f">{chromosme}\n") + pos + 1
        nucleotid = "T"

        index_result = self.reader.get_nucleotid(pos, chromosme)
        self.reader.fasta_file.seek(index_result, 0)
        nucleotid_result = self.reader.fasta_file.read(1)

        self.assertEqual(index, index_result)
        self.assertEqual(nucleotid, nucleotid_result)

    def test_get_nucleotid_last_line_prev_last(self):
        chromosme = "chr1"
        pos = 50 * 9 + 48
        index = len(f">{chromosme}\n") + pos + 9
        nucleotid = "A"

        index_result = self.reader.get_nucleotid(pos, chromosme)
        self.reader.fasta_file.seek(index_result, 0)
        nucleotid_result = self.reader.fasta_file.read(1)

        self.assertEqual(index, index_result)
        self.assertEqual(nucleotid, nucleotid_result)

    def test_get_nucleotid_last_line_last(self):
        chromosme = "chr1"
        pos = 50 * 9 + 49
        index = len(f">{chromosme}\n") + pos + 9
        nucleotid = "T"

        index_result = self.reader.get_nucleotid(pos, chromosme)
        self.reader.fasta_file.seek(index_result, 0)
        nucleotid_result = self.reader.fasta_file.read(1)

        self.assertEqual(index, index_result)
        self.assertEqual(nucleotid, nucleotid_result)

    def test_get_nucleotid_last_chr_first_line_first(self):
        chromosme = "chr3"
        pos = 0
        prev_elements = 10 * 51 * 2 + 6 * 2
        index = len(f">{chromosme}\n") + prev_elements + pos
        nucleotid = "T"

        index_result = self.reader.get_nucleotid(pos, chromosme)
        self.reader.fasta_file.seek(index_result, 0)
        nucleotid_result = self.reader.fasta_file.read(1)

        self.assertEqual(index, index_result)
        self.assertEqual(nucleotid, nucleotid_result)

    def test_get_nucleotid_last_chr_first_line_second(self):
        chromosme = "chr3"
        pos = 1
        prev_elements = 10 * 51 * 2 + 6 * 2
        index = len(f">{chromosme}\n") + prev_elements + pos
        nucleotid = "A"

        index_result = self.reader.get_nucleotid(pos, chromosme)
        self.reader.fasta_file.seek(index_result, 0)
        nucleotid_result = self.reader.fasta_file.read(1)

        self.assertEqual(index, index_result)
        self.assertEqual(nucleotid, nucleotid_result)

    def test_get_nucleotid_last_chr_last_line_prev_last(self):
        chromosme = "chr3"
        pos = 50 * 9 + 48
        prev_elements = 10 * 51 * 2 + 6 * 2 + 9
        index = len(f">{chromosme}\n") + prev_elements + pos
        nucleotid = "A"

        index_result = self.reader.get_nucleotid(pos, chromosme)
        self.reader.fasta_file.seek(index_result, 0)
        nucleotid_result = self.reader.fasta_file.read(1)

        self.assertEqual(index, index_result)
        self.assertEqual(nucleotid, nucleotid_result)

    def test_get_nucleotid_last_chr_last_line_last(self):
        chromosme = "chr3"
        pos = 50 * 9 + 49
        prev_elements = 10 * 51 * 2 + 6 * 2 + 9
        index = len(f">{chromosme}\n") + prev_elements + pos
        nucleotid = "T"

        index_result = self.reader.get_nucleotid(pos, chromosme)
        self.reader.fasta_file.seek(index_result, 0)
        nucleotid_result = self.reader.fasta_file.read(1)

        self.assertEqual(index, index_result)
        self.assertEqual(nucleotid, nucleotid_result)

    def test_get_from_interval_first_chr_first_line_first_5_elements(self):
        chromosme = "chr1"
        pos = 0
        length = 5
        index = self.reader.get_nucleotid(pos, chromosme)

        result = self.reader.get_from_interval(index, length)

        self.assertEqual(result, "TATCA")

    def test_get_from_interval_first_chr_first_line_last_5_elements(self):
        chromosme = "chr1"
        pos = 49 - 5 + 1
        length = 5
        index = self.reader.get_nucleotid(pos, chromosme)

        result = self.reader.get_from_interval(index, length)

        self.assertEqual(result, "AGGGC")

    def test_get_from_interval_first_chr_fourth_line_20_to_25_elements(self):
        chromosme = "chr1"
        pos = 50 * 3 + 19
        length = 5
        index = self.reader.get_nucleotid(pos, chromosme)

        result = self.reader.get_from_interval(index, length)

        self.assertEqual(result, "CTTCC")

    def test_get_from_interval_first_chr_fourth_line_48_to_next_2_elements(self):
        chromosme = "chr1"
        pos = 50 * 3 + 47
        length = 5
        index = self.reader.get_nucleotid(pos, chromosme)

        result = self.reader.get_from_interval(index, length)

        self.assertEqual(result, "GATGG")

    def test_get_from_interval_first_chr_last_line_first_5_elements(self):
        chromosme = "chr1"
        pos = 50 * 9
        length = 5
        index = self.reader.get_nucleotid(pos, chromosme)

        result = self.reader.get_from_interval(index, length)

        self.assertEqual(result, "CTGAT")

    def test_get_from_interval_first_chr_last_line_last_5_elements(self):
        chromosme = "chr1"
        pos = 50 * 9 + 49 - 5 + 1
        length = 5
        index = self.reader.get_nucleotid(pos, chromosme)

        result = self.reader.get_from_interval(index, length)

        self.assertEqual(result, "CTAAT")

    def test_get_prefix_first_row_first_element(self):
        chromosme = "chr1"
        pos = 0
        length = 5

        result = self.reader.get_prefix(pos, length, chromosme)

        self.assertEqual(result, "")

    def test_get_prefix_first_row_second_element(self):
        chromosme = "chr1"
        pos = 1
        length = 5

        result = self.reader.get_prefix(pos, length, chromosme)

        self.assertEqual(result, "T")

    def test_get_prefix_first_row_fifth_element(self):
        chromosme = "chr1"
        pos = 4
        length = 5

        result = self.reader.get_prefix(pos, length, chromosme)

        self.assertEqual(result, "TATC")

    def test_get_prefix_first_row_sixth_element(self):
        chromosme = "chr1"
        pos = 5
        length = 5

        result = self.reader.get_prefix(pos, length, chromosme)

        self.assertEqual(result, "TATCA")

    def test_get_prefix_second_row_first_element(self):
        chromosme = "chr1"
        pos = 50
        length = 5

        result = self.reader.get_prefix(pos, length, chromosme)

        self.assertEqual(result, "AGGGC")

    def test_get_prefix_second_row_second_element(self):
        chromosme = "chr1"
        pos = 51
        length = 5

        result = self.reader.get_prefix(pos, length, chromosme)

        self.assertEqual(result, "GGGCA")

    def test_get_prefix_first_element_next_chromosme(self):
        chromosme = "chr2"
        pos = 0
        length = 5

        result = self.reader.get_prefix(pos, length, chromosme)

        self.assertEqual(result, "")

    def test_get_prefix_more_that_one_line_chromosme(self):
        chromosme = "chr2"
        pos = 50 * 2 + 10
        first_line_sequence = "GTCGGAGGGC"
        second_line_sequence = "ATAACATCGTCGTGCTCCACGCTGTAATATCGTGCCAGGAAGTGACCGTG"
        thirsd_line_sequence = "AAGCAGCGTC"
        length = 10 + 50 + 10
        sequence = f"{first_line_sequence}{second_line_sequence}{thirsd_line_sequence}"

        result = self.reader.get_prefix(pos, length, chromosme)

        self.assertEqual(result, sequence)

    def test_get_sequence_sufix_first_line_first_element_5_length(self):
        chromosme = "chr2"
        pos = 0
        length = 5
        sequence = "ATCAA"

        result = self.reader.get_suffix(pos, length, chromosme)

        self.assertEqual(result, sequence)

    def test_get_sequence_sufix_first_line_prev_last_element_5_length(self):
        chromosme = "chr2"
        pos = 48
        length = 5
        sequence = "CATAA"

        result = self.reader.get_suffix(pos, length, chromosme)

        self.assertEqual(result, sequence)

    def test_get_sequence_sufix_more_than_one_line(self):
        chromosme = "chr2"
        pos = 52
        length = 47 + 50 + 10
        first_sequence = "ACATCGTCGTGCTCCACGCTGTAATATCGTGCCAGGAAGTGACCGTG"
        second_seqeunce = "AAGCAGCGTCGACTACTAGTTGATATCATCTATAGAACGCACGTGTATCT"
        last_sequence = "TTAAAGCGAA"
        sequence = f"{first_sequence}{second_seqeunce}{last_sequence}"

        result = self.reader.get_suffix(pos, length, chromosme)

        self.assertEqual(result, sequence)

    def test_get_nucleotide_sequence_firs_line_first_element(self):
        chromosme = "chr1"
        pos = 0

        result = self.reader.get_nucleotide_sequence(chromosme, pos)

        self.assertEqual(result, "T")

    def test_get_nucleotide_sequence_firs_line_lasts_6_elements(self):
        chromosme = "chr1"
        pos = 49
        length = 6

        result = self.reader.get_nucleotide_sequence(chromosme, pos, length=length)

        self.assertEqual(result, "CATAAC")

    def test_get_sequence(self):
        chromosome = "chr1"
        nucleotide = "TCGG"
        pos = 41
        from_nuc = 10
        prefix = "GCGAACCGCG"
        to_nuc = 20
        suffix = "AGGGCATAACATCGTCGTGC"

        result = self.reader.get_sequence(
            chromosome,
            nucleotide,
            pos,
            from_nuc,
            to_nuc,
        )

        self.assertEqual(result[0], prefix)
        self.assertEqual(result[1], nucleotide)
        self.assertEqual(result[2], suffix)
