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

    def generate_model_data(self, path: str):
        """Generates a file with the sequences mutated like:
            sequence:               ACGTGGTCAAGTCC
            sequence_simplified:    nnnnnnnmmmnnnn
        
        Parameters
        ----------
        path : str
            Path to store the data
        """

        #TODO: acabar método
        pass
