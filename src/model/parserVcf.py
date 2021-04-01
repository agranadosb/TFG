# -*- coding: utf-8 -*-

from src.vcf.vcfReader import VcfMutationsReader


class ParserVcf(object):
    """Parses data from vcf and fasta and prepare that data for a
    machine learning model

    Parameters
    ----------
    vcf_path : str
        Path of the vcf file
    fasta_path : str
        Path of the fasta file
    """

    def __init__(self, vcf_path: str, fasta_path: str):
        self.vcf_reader = VcfMutationsReader(vcf_path, fasta_path)

    def get_vcf(self):
        """Returns vcf file

        Returns
        -------
        str
            vcf file
        """
        return self.vcf_reader.get_vcf()

    def get_vcf_reader(self):
        """Returns vcf reader

        Returns
        -------
        VcfMutationsReader
            vcf reader
        """
        return self.vcf_reader

    def get_simplified_sequence(self, sequence: tuple):
        """Generates the simplified sequence by a equence:
            sequence:               (ACGTGGT,CAA,GTCC)
            sequence_simplified:    (lllllll,mmm,rrrr)

        Parameters
        ----------
        sequence : tuple
            Sequence to simplify

        Returns
        -------
        tuple
            sequence simplified
        """

        left_nucleotides = 'l' * len(sequence[0])
        mutations_nucleotides = 'm' * len(sequence[1])
        right_nucleotides = 'r' * len(sequence[2])

        return (left_nucleotides, mutations_nucleotides, right_nucleotides)

    def generate_simplified_sequences(self, path: str, filename: str = False, write_chromosme: bool = False):
        """Generates a file with the sequences mutated like:
            sequence:               ACGTGGTCAAGTCC
            sequence_simplified:    nnnnnnnmmmnnnn

        Parameters
        ----------
        path : str
            Path to store the data
        path : str, optional
            Filename of the result file
        write_chromosme : bool, optional
            Indicates the chromosme where the sequences being
        """

        if not filename:
            filename = 'parsed_data.pvcf'

        current_chromosme = ''
        with open(f'{path}/{filename}', 'w') as parsed_data_file:
            for i in self.get_vcf():
                seq = self.vcf_reader.get_sequence(i.CHROM, i.REF, i.POS, 5, 5)

                sequence_label = "".join(seq)
                simplifed_sequence_label = "".join(self.get_simplified_sequence(seq))

                prefix = ''
                if write_chromosme:
                    prefix = f'{i.CHROM}\t'

                parsed_data_file.write(
                    f'{prefix}{sequence_label}\n{prefix}{simplifed_sequence_label}\n'
                )
