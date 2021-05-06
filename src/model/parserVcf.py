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
        return (sequence[0], mutation[1].sequence.lower(), sequence[2])

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
            "l" * len(sequence[0]),
            "m" * len(mutation[1].sequence),
            "r" * len(sequence[2]),
        )

    def get_extended_sequence(self, sequence: tuple, mutation: str):
        """Generates the simplified sequence by a equence:
            sequence:               (ACGTGGT,CAA,GTCC)
            sequence_simplified:    ([1,2,3,4,3,3,4],[12,11,11],[23,24,22,22])

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

        left = [self.left_map(nucletid) for nucletid in sequence[0]]
        middle = [self.middle_map(nucletid) for nucletid in sequence[1].sequence]
        right = [self.right_map(nucletid) for nucletid in sequence[2]]

        return (left, middle, right)

    def generate_sequences(
        self,
        path: str,
        filename: str = False,
        write_chromosme: bool = False,
        method=False,
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
        method : function, optional
            Method to generate the mutated sequences, simplified as default
        add_original : bool, optional
            If true adds the original sequence into the file
        """

        if not method:
            method = self.get_simplified_sequence

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
                    f'{original_sequence}{prefix}{"".join(method(sequence, i.ALT))}\n'
                )

    def generate_lower_sequences(
        self,
        path: str,
        filename: str = False,
        write_chromosme: bool = False,
        add_original: bool = True,
    ):
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

        self.generate_sequences(
            path,
            filename or "parsed_lower_data.pvcf",
            write_chromosme,
            method=self.get_lower_sequence,
            add_original=add_original,
        )

    def generate_simplified_sequences(
        self,
        path: str,
        filename: str = False,
        write_chromosme: bool = False,
        add_original: bool = True,
    ):
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

        self.generate_sequences(
            path,
            filename or "parsed_simplified_data.pvcf",
            write_chromosme,
            add_original=add_original,
        )

    def generate_extended_sequences(
        self,
        path: str,
        filename: str = False,
        write_chromosme: bool = False,
        add_original: bool = True,
    ):
        """Generates a file with the sequences mutated like:
            sequence:               A-C-G-T-G-G-T- C- A- A- G- T- C- C
            sequence_simplified:    1-2-3-4-3-3-4-12-11-11-23-24-22-22

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

        letters = ["A", "C", "G", "T"]

        self.left_map = {key: value for (key, value) in zip(letters, range(4))}
        self.middle_map = {key: value + 10 for (key, value) in zip(letters, range(4))}
        self.right_map = {key: value + 20 for (key, value) in zip(letters, range(4))}

        self.generate_sequences(
            path,
            filename or "parsed_extended_data.pvcf",
            write_chromosme,
            method=self.get_extended_sequence,
            add_original=add_original,
        )

    """ TODO: Separar en clases según el tipo de parser, y heredar la clase común con los métodos
        que deberá ser abstract
    """
