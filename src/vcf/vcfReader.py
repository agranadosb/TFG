# -*- coding: utf-8 -*-

from src.lib.PyVCF.vcf import Reader as VcfReader
import gzip
import os
import shutil
import sys


class VcfMutationsReader(object):

    def __init__(self, vcf_path: str, fasta_path: str):
        self.vcf_file = VcfReader(open(vcf_path, 'r'))
        self.fasta_file = gzip.open(fasta_path, 'r')
        self.fatsa_indexes = {}
        self.fasta_filename = fasta_path.replace('.gz', '')
        self.fasta_filename_index = fasta_path.replace('.gz', '.index')
        self.fasta_file_line_length = 50

        with open(self.fasta_filename, 'wb') as f_out:
            shutil.copyfileobj(self.fasta_file, f_out)

        self.fasta_file = open(self.fasta_filename, 'r')

        self.generateFastaIndex()

        for i in self.fasta_file:
            if not i.startswith('>'):
                self.fasta_file_line_length = len(i) - 1
                break

    def getVcf(self):
        return self.vcf_file

    def getFasta(self):
        return self.fasta_file

    def getFastaIndexes(self):
        return self.fatsa_indexes

    def generateFastaIndex(self):
        os.system(
            f"cat {self.fasta_filename} | grep -n '>' > {self.fasta_filename_index}")
        with open(self.fasta_filename_index, 'r') as fasta_index_file:
            for i in fasta_index_file:
                index_id = i.split(':')

                self.fatsa_indexes[
                    index_id[1].rstrip().replace('>', '')
                ] = int(index_id[0]) - 1
        self.fasta_keys = list(self.fatsa_indexes.keys())

    def _get_seq_interval(self, from_seq_nuc: int, from_index: int):
        from_seq = ''
        if from_seq_nuc:
            self.fasta_file.seek(from_index, 0)

            rest_from = from_seq_nuc
            while rest_from > 0:
                read_seq = self.fasta_file.readline().rstrip()[:rest_from]
                from_seq += read_seq
                rest_from -= len(read_seq)

        return from_seq

    def get_cromosme_index(self, cromosme_number: str):
        cromosme_index = 0
        if self.fatsa_indexes[cromosme_number] > 0:
            cromosme_index = self.fasta_keys.index(cromosme_number)
            cyhrs_string_sum = sum(
                [len(self.fasta_keys[i]) + 2 for i in range(cromosme_index)])

            lines = self.fatsa_indexes[cromosme_number] - cromosme_index
            lines_sum = lines * (self.fasta_file_line_length + 1)

            cromosme_index = lines_sum + cyhrs_string_sum

        return cromosme_index

    def _get_seq_prefix(self, line_start, pointer_index, from_nuc, last_line_nuc, line_index_ends, cromosme_number_length):
        from_nuc_end = pointer_index - from_nuc - 1 + last_line_nuc
        if line_index_ends - pointer_index <= self.fasta_file_line_length and not last_line_nuc:
            from_nuc_end += 1
        if from_nuc_end - (line_start + cromosme_number_length) <= cromosme_number_length:
            from_nuc = from_nuc_end - (line_start + cromosme_number_length)
            if from_nuc <= 0:
                from_nuc += cromosme_number_length
                from_nuc_end = line_start + cromosme_number_length + 1
        return self._get_seq_interval(from_nuc, from_nuc_end)

    def _get_seq_sufix(self, pointer_index, to_nuc, line_index_ends, next_cromosome_length, file_ends):
        to_nuc_end = pointer_index + 1

        to_nuc_cpoy = to_nuc
        if to_nuc_end + to_nuc >= line_index_ends:
            to_nuc = line_index_ends - (to_nuc_end + to_nuc) + \
                next_cromosome_length - 1
            if to_nuc <= 0 and file_ends:
                to_nuc += to_nuc_cpoy

        return self._get_seq_interval(to_nuc, pointer_index + 1)

    def _get_prefix_and_sufix(self, cromosme_number, line_start, pointer_index, from_nuc, to_nuc, last_line_nuc, cromosme_number_length):
        next_cromosome_index = self.fasta_keys.index(cromosme_number) + 1
        line_index_ends = False
        next_cromosome_length = False
        file_ends = False

        if next_cromosome_index == len(self.fasta_keys):
            line_index_ends = os.stat(self.fasta_filename).st_size
            next_cromosome_length = 0
            file_ends = True
        elif next_cromosome_index < len(self.fasta_keys):
            line_index_ends = self.get_cromosme_index(
                self.fasta_keys[next_cromosome_index])
            next_cromosome_length = len(
                f'>{self.fasta_keys[next_cromosome_index]}')

        from_seq = self._get_seq_prefix(
            line_start,
            pointer_index,
            from_nuc,
            last_line_nuc,
            line_index_ends,
            cromosme_number_length
        )

        to_seq = self._get_seq_sufix(
            pointer_index,
            to_nuc,
            line_index_ends,
            next_cromosome_length,
            file_ends
        )

        return {
            'prefix': from_seq,
            'sufix': to_seq,
        }

    def get_seq_by_chr_pos(self, cromosme_number: str, pos: int, from_nuc: int = False, to_nuc: int = False,):
        # Get the start of the chromosme in the fasta file
        line_start = self.get_cromosme_index(cromosme_number)
        cromosme_number_length = len(f'>{cromosme_number}')
        num_new_lines = int(pos / self.fasta_file_line_length)
        last_line_nuc = pos % self.fasta_file_line_length == 0

        # Get the position of the nucleotid on the file (1 char is one byte, that's why we use seek)
        pointer_index = pos + line_start + cromosme_number_length + num_new_lines - last_line_nuc

        self.fasta_file.seek(pointer_index, 0)

        nucleotide = self.fasta_file.readline()[0]

        if not from_nuc and not to_nuc:
            return nucleotide

        prefix_and_sufix = self._get_prefix_and_sufix(
            cromosme_number,
            line_start,
            pointer_index,
            from_nuc, to_nuc,
            last_line_nuc,
            cromosme_number_length
        )

        from_seq = prefix_and_sufix['prefix']
        to_seq = prefix_and_sufix['sufix']

        return (from_seq, nucleotide, to_seq)
