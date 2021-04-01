# -*- coding: utf-8 -*-

from src.lib.PyVCF.vcf import Reader as VcfReader
import gzip
import os
import shutil
import sys


class VcfMutationsReader(object):
    """Gets information of mutations by a vcf and a fasta files

    Using a fasta file and a vcf file gets information and parses data
    from fasta file using the mutations on vcf

    Parameters
    ----------
    vcf_path : str
        Path of the vcf file
    fasta_path : str
        Path of the fasta file
    """

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

        self.generate_fasta_information()

        for i in self.fasta_file:
            if not i.startswith('>'):
                self.fasta_file_line_length = len(i) - 1
                break

    def get_vcf(self):
        """Returns vcf file

        Returns
        -------
        str
            vcf file
        """
        return self.vcf_file

    def get_fasta(self):
        """Returns fasta file

        Returns
        -------
        str
            fasta file
        """
        return self.fasta_file

    def get_fatsa_chromosme_information(self):
        """Returns information of chromosomes of fasta file

        Returns
        -------
        dict
            fasta information file
        """
        return self.fatsa_chromosme_information

    def set_chromosme_information(self, chromosme):
        """Adds chromosome information of 'chromosme'

        Adds chromosome information of 'chromosme':
         - label_length: length of the chromosome label
         - file_line: the line where the chromosome starts in fasta file
         - line_start: index where chromosome is in fasta file

        Parameters
        ----------
        chromosme : str
            Chromosome label as '>chromosme:line_fasta'
        """
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
        """Adds chromosome information about sizes in fasta file

        Adds chromosome information about sizes in fasta file:
         - line_index_ends: where the chromosome ends
         - next_cromosome_length: label length of the next chromosome
         - file_ends: shows if the chromosome is the last
         - index: chromosome

        Parameters
        ----------
        chromosme : str
            Chromosome label as '>chromosme:line_fasta'
        """
        chromosome_information = self.fatsa_chromosme_information[self.fasta_keys[chromosme_index]]
        line_index_ends = os.stat(self.fasta_filename).st_size
        next_cromosome_length = 0
        file_ends = True

        if chromosme_index != len(self.fasta_keys) - 1:
            next_chromosome = self.fasta_keys[chromosme_index + 1]
            line_index_ends = self.fatsa_chromosme_information[next_chromosome]['line_start']
            next_cromosome_length = len(f'>{next_chromosome}')
            file_ends = False

        chromosome_information['line_index_ends'] = line_index_ends
        chromosome_information['next_cromosome_length'] = next_cromosome_length
        chromosome_information['file_ends'] = file_ends
        chromosome_information['index'] = chromosme_index

    def generate_fasta_information(self):
        """Generates fasta file information about the chromosomes

        Generates fasta file information about the chromosomes:
         - label_length: length of the chromosome label
         - file_line: the line where the chromosome starts in fasta file
         - line_start: index where chromosome is in fasta file
         - line_index_ends: where the chromosome ends
         - next_cromosome_length: label length of the next chromosome
         - file_ends: shows if the chromosome is the last
         - index: chromosome
        """
        command = f"cat {self.fasta_filename} | grep -n '>' > {self.fasta_filename_index}"
        os.system(command)
        with open(self.fasta_filename_index, 'r') as fasta_index_file:
            for i in fasta_index_file:
                self.set_chromosme_information(i)

        for i in range(len(self.fasta_keys)):
            self.set_chromosme_sizes(i)

    def get_seq_interval(self, from_seq_nuc: int, from_index: int):
        """Returns the sequence that starts in 'from_index' and has 'from_seq_nuc' length

        Parameters
        ----------
        from_seq_nuc : int
            Number of elements from the sequence
        from_index : int
            Start of the sequence

        Returns
        -------
        str
            sequence

        Raises
        ------
        IndexError
            When the sequence from the 'from_index' is over the
            end of the file or is before the start of the file
        """
        from_seq = ''
        prev = -1
        if from_seq_nuc:
            self.fasta_file.seek(from_index, 0)

            rest_from = from_seq_nuc
            while rest_from > 0:
                read_seq = self.fasta_file.readline().rstrip()[:rest_from]
                from_seq += read_seq
                rest_from -= len(read_seq)
                if prev == rest_from:
                    raise IndexError('wrong indexes')
                prev = rest_from

        return from_seq

    def get_chromosome_file_index(self, chromosome_label: str):
        """Returns index of the chromosome 'chromosome_label' on the fasta file

        Parameters
        ----------
        chromosome_label : int
            Label of the chromosome

        Returns
        -------
        str
            chromosome index on fasta file
        """
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

    def get_seq_prefix(self, pos, from_nuc, last_line_nuc, chromosome_label):
        """Gets prefix from an index on the fasta file

        Parameters
        ----------
        pos : int
            Index on chromosome 'chromosome_label'
        from_nuc : int
            Length of the prefix
        last_line_nuc : bool
            Indicates if the index is in the last line of a chromosome
        chromosome_label : str
            Label of the chromosome

        Returns
        -------
        str
            prefix
        """
        pointer_index = self.get_nucleotid_fasta_index(pos, chromosome_label)
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

    def get_seq_sufix(self, pos, to_nuc, chromosome_label):
        """Gets suffix from an index on the fasta file

        Parameters
        ----------
        pos : int
            Index on chromosome 'chromosome_label'
        to_nuc : int
            Length of the suffix
        chromosome_label : str
            Label of the chromosome

        Returns
        -------
        str
            suffix
        """
        pointer_index = self.get_nucleotid_fasta_index(pos, chromosome_label)
        file_ends = self.fatsa_chromosme_information[chromosome_label]['file_ends']
        next_cromosome_length = self.fatsa_chromosme_information[
            chromosome_label
        ]['next_cromosome_length']
        line_index_ends = self.fatsa_chromosme_information[chromosome_label]['line_index_ends']

        to_nuc_end = pointer_index + 1
        to_nuc_cpoy = to_nuc

        if to_nuc_end + to_nuc >= line_index_ends:
            to_nuc = line_index_ends - (to_nuc_end + to_nuc) + \
                next_cromosome_length - 1
            if to_nuc <= 0 and file_ends:
                to_nuc += to_nuc_cpoy

        return self.get_seq_interval(to_nuc, pointer_index + 1)

    def get_nucleotid_fasta_index(self, pos, chromosome_label):
        """Gets the index of a nucleotide by its position in a chromosome

        Parameters
        ----------
        pos : int
            Index on chromosome 'chromosome_label'
        chromosome_label : str
            Label of the chromosome

        Returns
        -------
        str
            index of the nucleotid in the fasta file
        """
        # It's necessary to taking into account that seek method counts new lines as a character
        num_new_lines = int(pos / self.fasta_file_line_length)
        last_line_nuc = pos % self.fasta_file_line_length == 0
        line_start = self.fatsa_chromosme_information[chromosome_label]['line_start']
        label_length = self.fatsa_chromosme_information[chromosome_label]['label_length']

        # Get the position of the nucleotid on the file (1 char is one byte, that's why we use seek)
        return pos + line_start + label_length + num_new_lines - last_line_nuc

    def get_nucleotide(self, chromosome_label: str, pos: int, length: int = 1):
        """Gets the nucleotide (or sequence) from the psoition 'pos' of a chromosme

        Parameters
        ----------
        chromosome_label : str
            The chromosome where the nucleotide is going to be obtained
        pos : int
            Position of the nucleotide on the chromosome 'chromosome_label'
        length : str, optional
            If its a sequence of nucleotides what it wants to be obtained,
            indicates the length of that sequence

        Returns
        -------
        str
            nucleotid in the fasta
        """

        nucleotide = ''
        pointer_index = self.get_nucleotid_fasta_index(pos, chromosome_label)
        for i in range(length):
            self.fasta_file.seek(pointer_index, 0)
            readed_character = self.fasta_file.read(1)

            if readed_character == '\n':
                readed_character = self.fasta_file.read(1)

            nucleotide += readed_character
            pointer_index += 1

        return nucleotide

    def get_sequence(self, chromosome_label: str, nucleotide_label: str, pos: int, from_nuc: int, to_nuc: int,):
        """Gets the a fasta's sequence by a position in a chromosome.

        Gets the nucleotide in position 'pos' from the chromosme 'chromosome_label'
        and, if the parameters 'from_nuc' or 'to_nuc' are given, the prefix and suffix 
        with 'from_nuc' and 'to_nuc' length respectively as tuple (prefix, nucleotid, suffix)

        Parameters
        ----------
        chromosome_label : str
            The chromosome where the nucleotide is going to be obtained
        nucleotide_label : str
            Label of the nucleotide
        pos : int
            Position of the nucleotide on the chromosome 'chromosome_label'
        from_nuc : int
            Length of the prefix of the sequence
        to_nuc : int
            Length of the suffix of the sequence

        Returns
        -------
        str
            nucleotid in the fasta
        tuple
            nucleotide with the prefix and the suffix
        """
        last_line = pos % self.fasta_file_line_length == 0
        nucleotide_length = len(nucleotide_label)

        pref = self.get_seq_prefix(pos, from_nuc, last_line, chromosome_label)
        nucleotide = self.get_nucleotide(chromosome_label, pos, nucleotide_length)
        suff = self.get_seq_sufix(pos + nucleotide_length - 1, to_nuc, chromosome_label)

        return (pref, nucleotide, suff)
