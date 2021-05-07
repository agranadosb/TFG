# -*- coding: utf-8 -*-

from abc import ABC, abstractmethod

from src.vcf.vcfReader import VcfMutationsReader


class ParserVcf(ABC):
    """Parses data from vcf and fasta and prepare that data for a
    machine learning model

    Parameters
    ----------
    vcf_path : str
        Path of the vcf file
    fasta_path : str
        Path of the fasta file
    """

    name = False

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

    @abstractmethod
    def method(self, sequence: tuple, mutation: str):
        pass

    def generate_sequences(
        self,
        path: str,
        filename: str = False,
        write_chromosme: bool = False,
        add_original: bool = True,
    ):
        """Generates a file with the sequences mutated using 'method'

        Parameters
        ----------
        path : str
            Path to store the data
        path : str, optional
            Filename of the result file
        write_chromosme : bool, optional
            Indicates the chromosme where the sequences being
        add_original : bool, optional
            If true adds the original sequence into the file
        """

        if not filename:
            filename = f"parsed_{self.name}_data.pvcf"

        with open(f"{path}/{filename}", "w") as parsed_data_file:
            for i in self.get_vcf():
                sequence = self.vcf_reader.get_sequence(i.CHROM, i.REF, i.POS, 5, 5)

                prefix = ""
                if write_chromosme:
                    prefix = f"{i.CHROM}\t"

                original_sequence = ""
                if add_original:
                    original_sequence = (
                        f'{prefix}{"".join(sequence)}\n' if add_original else ""
                    )

                parsed_data_file.write(
                    f'{original_sequence}{prefix}{"".join(self.method(sequence, i.ALT))}\n'
                )
