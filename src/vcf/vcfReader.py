# -*- coding: utf-8 -*-

from src.lib.PyVCF.vcf import Reader as VcfReader
import gzip
import os
import shutil


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
            chrs = list(self.fatsa_indexes.keys())
            cyhrs_string_sum = sum(
                [len(chrs[i]) + 2 for i in range(chrs.index(cromosme_number))])
            cromosme_index = (
                self.fatsa_indexes[cromosme_number] - 1) * (self.fasta_file_line_length + 1) + cyhrs_string_sum

        return cromosme_index

    def get_seq_by_chr_pos(self, cromosme_number: str, pos: int, from_nuc: int = False, to_nuc: int = False,):
        # Get the start of the chromosme in the fasta file
        line_start = self.get_cromosme_index(cromosme_number)

        pointer_index = pos + line_start + len(f'>{cromosme_number}')

        num_new_lines = int(pos / self.fasta_file_line_length)
        last_line = pos % self.fasta_file_line_length == 0

        # Get the position of the nucleotid on the file (1 char is one byte, that's why we use seek)
        pointer_index += num_new_lines - last_line

        self.fasta_file.seek(pointer_index, 0)

        nucleotide = self.fasta_file.readline()[0]

        """ TODO:
            Hay que tener en cuenta que si se busca un nucleotido y su contexto (n por delante y m por detras),
            si n o m supera el final del cromosma hay que devolver sólo el principio o el final del cromosoma
            respectivamente

            >chrX
            ACGTAAGGT*C*CAGTTGCAAAA
            >chrY
            ...

            con n = 12 y m = 20 tiene que devolver:
                ('ACGTAAGGT', 'C' 'CAGTTGCAAAA')
            
            Es decir, cogiendo el max(lo que queda por delante, n) y max(lo que queda por detrás, m)

            Esto se puede hacer con el índice del cromosoma en el que se está buscando actualmente y el
            siguiente (obteniendo los índices en el archivo del cromosoma)
        """
        from_seq = self._get_seq_interval(
            from_nuc, pointer_index - from_nuc - 1 + last_line)
        to_seq = self._get_seq_interval(to_nuc, pointer_index + 1)

        if from_seq or to_seq:
            return (from_seq, nucleotide, to_seq)

        return nucleotide
