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
        self.fatsa_chromosme_information = {}
        self.fasta_keys = []
        self.fasta_filename = fasta_path.replace('.gz', '')
        self.fasta_filename_index = fasta_path.replace('.gz', '.index')
        self.fasta_file_line_length = 50

        with open(self.fasta_filename, 'wb') as f_out:
            shutil.copyfileobj(self.fasta_file, f_out)

        self.fasta_file = open(self.fasta_filename, 'r')

        self.generate_fasta_index()

        for i in self.fasta_file:
            if not i.startswith('>'):
                self.fasta_file_line_length = len(i) - 1
                break

    def get_vcf(self):
        return self.vcf_file

    def get_fasta(self):
        return self.fasta_file

    def get_fatsa_chromosme_information(self):
        return self.fatsa_chromosme_information

    def set_chromosme_information(self, chromosme):
        index_id = chromosme.split(':')
        chromosme_label = index_id[1].rstrip().replace('>', '')

        self.fasta_keys.append(chromosme_label)

        # Add chromosme information
        self.fatsa_chromosme_information[chromosme_label] = {
            'label_length': len(f'>{chromosme_label}'),
            'file_line': int(index_id[0]) - 1
        }

        self.fatsa_chromosme_information[
            chromosme_label
        ]['line_start'] = self.get_chromosome_file_index(chromosme_label)

    def set_chromosme_sizes(self, chromosme_index):
        chromosome_information = self.fatsa_chromosme_information[self.fasta_keys[chromosme_index]]
        line_index_ends = os.stat(self.fasta_filename).st_size
        next_cromosome_length = 0
        file_ends = True

        if chromosme_index != len(self.fasta_keys) - 1:
            next_chromosome = self.fasta_keys[chromosme_index + 1]
            line_index_ends = self.get_chromosome_file_index(next_chromosome)
            next_cromosome_length = len(f'>{next_chromosome}')
            file_ends = False

        chromosome_information['line_index_ends'] = line_index_ends
        chromosome_information['next_cromosome_length'] = next_cromosome_length
        chromosome_information['file_ends'] = file_ends
        chromosome_information['index'] = chromosme_index

    def generate_fasta_index(self):
        command = f"cat {self.fasta_filename} | grep -n '>' > {self.fasta_filename_index}"
        os.system(command)
        with open(self.fasta_filename_index, 'r') as fasta_index_file:
            for i in fasta_index_file:
                self.set_chromosme_information(i)

        for i in range(len(self.fasta_keys)):
            self.set_chromosme_sizes(i)

    def get_seq_interval(self, from_seq_nuc: int, from_index: int):
        from_seq = ''
        if from_seq_nuc:
            self.fasta_file.seek(from_index, 0)

            rest_from = from_seq_nuc
            while rest_from > 0:
                read_seq = self.fasta_file.readline().rstrip()[:rest_from]
                from_seq += read_seq
                rest_from -= len(read_seq)

        return from_seq

    def get_chromosome_file_index(self, chromosome_label: str):
        file_line = self.fatsa_chromosme_information[chromosome_label]['file_line']

        if file_line <= 0:
            return 0

        chromosome_index = self.fasta_keys.index(chromosome_label)
        cyhrs_string_sum = sum(
            [len(self.fasta_keys[i]) + 2 for i in range(chromosome_index)])

        lines_sum = (file_line - chromosome_index) * \
            (self.fasta_file_line_length + 1)

        chromosome_file_index = lines_sum + cyhrs_string_sum

        return chromosome_file_index

    def get_seq_prefix(self, pointer_index, from_nuc, last_line_nuc, chromosome_label):
        line_start = self.fatsa_chromosme_information[chromosome_label]['line_start']
        chromosome_label_length = self.fatsa_chromosme_information[chromosome_label]['label_length']
        line_index_ends = self.fatsa_chromosme_information[chromosome_label]['line_index_ends']

        from_nuc_end = pointer_index - from_nuc - 1 + last_line_nuc

        if line_index_ends - pointer_index <= self.fasta_file_line_length and not last_line_nuc:
            from_nuc_end += 1

        if from_nuc_end - (line_start + chromosome_label_length) <= chromosome_label_length:
            from_nuc = from_nuc_end - (line_start + chromosome_label_length)
            if from_nuc <= 0:
                from_nuc += chromosome_label_length
                from_nuc_end = line_start + chromosome_label_length + 1

        return self.get_seq_interval(from_nuc, from_nuc_end)

    def get_seq_sufix(self, pointer_index, to_nuc, chromosome_number):
        file_ends = self.fatsa_chromosme_information[chromosome_number]['file_ends']
        next_cromosome_length = self.fatsa_chromosme_information[
            chromosome_number
        ]['next_cromosome_length']
        line_index_ends = self.fatsa_chromosme_information[chromosome_number]['line_index_ends']

        to_nuc_end = pointer_index + 1
        to_nuc_cpoy = to_nuc

        if to_nuc_end + to_nuc >= line_index_ends:
            to_nuc = line_index_ends - (to_nuc_end + to_nuc) + \
                next_cromosome_length - 1
            if to_nuc <= 0 and file_ends:
                to_nuc += to_nuc_cpoy

        return self.get_seq_interval(to_nuc, pointer_index + 1)

    def get_nucleotid_fasta_index(self, pos, chromosome_number):
        # It's necessary to taking into account that seek method counts new lines as a character
        num_new_lines = int(pos / self.fasta_file_line_length)
        last_line_nuc = pos % self.fasta_file_line_length == 0
        line_start = self.fatsa_chromosme_information[chromosome_number]['line_start']
        label_length = self.fatsa_chromosme_information[chromosome_number]['label_length']

        # Get the position of the nucleotid on the file (1 char is one byte, that's why we use seek)
        pointer_index = pos + line_start + label_length + num_new_lines - last_line_nuc

        return pos + line_start + label_length + num_new_lines - last_line_nuc

    def get_seq_by_chr_pos(self, chromosome_number: str, pos: int, from_nuc: int = False, to_nuc: int = False,):
        """Gets the a fasta's sequence by a position in a chromosome.

        Gets the nucleotide in position 'pos' from the chromosme 'chromosome_number'
        and, if the parameters 'from_nuc' or 'to_nuc' are given, the prefix and suffix 
        with 'from_nuc' and 'to_nuc' length respectively as tuple (prefix, nucleotid, suffix)

        Parameters
        ----------
        chromosome_number : str
            The chromosome where the nucleotide is going to be obtained
        pos : int
            Position of the nucleotide on the chromosome 'chromosome_number'
        from_nuc : int, optional
            Length of the prefix of the sequence
        to_nuc : int, optional
            Length of the suffix of the sequence

        Returns
        -------
        str
            nucleotid in the fasta
        tuple
            nucleotide with the prefix and the suffix
        """

        pointer_index = self.get_nucleotid_fasta_index(pos, chromosome_number)

        self.fasta_file.seek(pointer_index, 0)

        nucleotide = self.fasta_file.readline()[0]

        if not from_nuc and not to_nuc:
            return nucleotide

        from_seq = self.get_seq_prefix(
            pointer_index,
            from_nuc,
            pos % self.fasta_file_line_length == 0,
            chromosome_number
        )

        to_seq = self.get_seq_sufix(
            pointer_index,
            to_nuc,
            chromosome_number
        )

        return (from_seq, nucleotide, to_seq)
