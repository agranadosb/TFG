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
        self.nucletodies_mutations_map = {
            'A': 'a',
            'C': 'c',
            'G': 'g',
            'T': 't',
        }

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

    def get_lower_sequence(self, sequence: tuple, mutation: str):
        """Generates the same sequence but in a lower case in the mutations part:
            sequence:               (ACGTGGT,CAA,GTCC)
            sequence_simplified:    (ACGTGGT,caa,GTCC)

        Parameters
        ----------
        sequence : tuple
            Sequence to simplify
        mutation : tuple
            Mutation sequence

        Returns
        -------
        tuple
            lower sequence
        """
        return (sequence[0], mutation[0].sequence.lower(), sequence[2])

    def get_simplified_sequence(self, sequence: tuple, mutation: str):
        """Generates the simplified sequence by a equence:
            sequence:               (ACGTGGT,CAA,GTCC)
            sequence_simplified:    (lllllll,mmm,rrrr)

        Parameters
        ----------
        sequence : tuple
            Sequence to simplify
        mutation : tuple
            Mutation sequence

        Returns
        -------
        tuple
            sequence simplified
        """

        return (
            'l' * len(sequence[0]),
            'm' * len(mutation),
            'r' * len(sequence[2])
        )

    def generate_sequences(self, path: str, filename: str = False, write_chromosme: bool = False, method=False, add_original: bool = True):
        """Generates a file with the sequences mutated using 'method'

        Parameters
        ----------
        path : str
            Path to store the data
        path : str, optional
            Filename of the result file
        write_chromosme : bool, optional
            Indicates the chromosme where the sequences being
        method : function, optional
            Method to generate the mutated sequences, simplified as default
        add_original : bool, optional
            If true adds the original sequence into the file
        """

        if not method:
            method = self.get_simplified_sequence

        current_chromosme = ''
        with open(f'{path}/{filename}', 'w') as parsed_data_file:
            for i in self.get_vcf():
                equence = self.vcf_reader.get_sequence(
                    i.CHROM,
                    i.REF,
                    i.POS,
                    5,
                    5
                )

                sequence_label = "".join(equence)
                parsed_sequence = "".join(method(equence, i.ALT))

                prefix = ''
                if write_chromosme:
                    prefix = f'{i.CHROM}\t'
                
                original_sequence = ''
                if add_original:
                    original_sequence = f'{prefix}{sequence_label}\n'

                mutated_sequence = f'{prefix}{parsed_sequence}\n'

                parsed_data_file.write(f'{original_sequence}{mutated_sequence}')

    def generate_lower_sequences(self, path: str, filename: str = False, write_chromosme: bool = False, add_original: bool = True):
        """Generates a file with the sequences mutated like:
            sequence:               ACGTGGTCAAGTCC
            sequence_simplified:    ACGTGGTcaaGTCC

        Parameters
        ----------
        path : str
            Path to store the data
        path : str, optional
            Filename of the result file
        write_chromosme : bool, optional
            Indicates the chromosme where the sequences being
        add_original: bool, optional
            If true adds the original sequence into the file
        """

        if not filename:
            filename = 'parsed_lower_data.pvcf'

        self.generate_sequences(
            path,
            filename,
            write_chromosme,
            method=self.get_lower_sequence,
            add_original=add_original
        )

    def generate_simplified_sequences(self, path: str, filename: str = False, write_chromosme: bool = False, add_original: bool = True):
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
        add_original: bool, optional
            If true adds the original sequence into the file
        """

        if not filename:
            filename = 'parsed_simplified_data.pvcf'

        self.generate_sequences(
            path,
            filename,
            write_chromosme,
            add_original=add_original
        )
