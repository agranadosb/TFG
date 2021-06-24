# -*- coding: utf-8 -*-

import gzip
import logging
import os
import shutil

from src.fasta.chromosome import Chromosome
from src.logging.tqdmLoggingHandler import TqdmLoggingHandler
from tqdm import tqdm


class FastaReader(object):
    """Gets a FASTA file and parse it, getting information about the chromosomes.

    After parse the fasta file, the class will have a dictionary with information per
    chromosome. Each value is a `src.fasta.chromosome.Chromosome` class.

    For example, after analize a Fasta file with two chromosomes, the dictionary will be:

    ```python
        chromosomes = {
            'chr1': Chromosome(...),
            'chr2': Chromosome(...)
        }
    ```

    This class allow to get sequences from the fasta file from chromosomes using
    `get_sequence` method.

    To acces to the information of one chromosome we can use array acces way:

    ```python
        fr = FastaReader(...)

        chr = fr['chr1']

        sequence = chr.sequence(...)
        sequence_string = chr[2:9]
    ```

    Parameters
    ----------
    fasta_path : str
        Path of the fasta file.
    """

    chromosomes: dict = {}
    """The dictionary that contains all chromosomes data of the FASTA file, Each key is
    the identifier of each chromosome."""

    _labels: dict = {}
    """Mapping between labels and chromosome identifier."""

    _chromosomes_list: list = []
    """List of the identifiers of all the chromosomes of the FASTA file."""

    fasta_filename: str = None
    """Name of the fasta file"""

    def __init__(self, fasta_path: str):
        logging.info("Loading fasta file")
        self.fasta_filename = fasta_path.replace(".gz", "")

        logging.info("Loading fasta information")
        self._set_fasta_file(fasta_path)
        self._parse_chromosomes()

    def _set_fasta_file(self, fasta_path: str):
        """Gets the fasta file from gz file and unzip the file to get the complete fasta
        file. If there is a unziped file with the same name without .gz extension, the
        class use that file instead of unzip the gz one.

        Parameters
        ----------
        fasta_path : str
            Path of the fasta file.
        """
        if not os.path.isfile(self.fasta_filename):
            logging.info("Unzip fasta gz file")
            self.fasta_file = gzip.open(fasta_path, "r")

            with open(self.fasta_filename, "wb") as f_out:
                shutil.copyfileobj(self.fasta_file, f_out)

        self.fasta_file = open(self.fasta_filename, "r")

    def _parse_chromosomes(self) -> dict:
        """Create all Cromosmes of he FASTA file and append it to a dictionary with the
        names of chromosomes as keys.

        For example, after analize a Fasta file with two chromosomes, the dictionary
        will be:

        ```python
            chromosomes = {
                'chr1': Chromosome(),
                'chr2': Chromosome()
            }
        ```

        Returns
        -------
        Chromosomes.
        """
        self._chromosomes_list = []
        self.chromosomes = {}
        logger = logging.getLogger()
        tqdm_out = TqdmLoggingHandler(logger, level=logging.INFO)

        lines = 0
        index = 0
        length = 0
        result = {}
        current_chromosome = False
        self.fasta_file.seek(0, 0)
        for i in tqdm(self.fasta_file, file=tqdm_out):
            lines += 1
            if i.startswith(">"):
                if current_chromosome:
                    result[current_chromosome]["length"] = length
                length = 0
                label_length = len(i)
                labels = i.rstrip().replace(">", "").split()
                current_chromosome = labels[len(labels) - 1]

                self._chromosomes_list.append(current_chromosome)
                result[current_chromosome] = {
                    "name": current_chromosome,
                    "line_length": False,
                    "label_length": label_length,
                    "index_start": index,
                    "length": False,
                    "labels": labels,
                }
                for label in labels:
                    self._labels[label] = current_chromosome
                is_populated = False
                index += len(i)
                continue

            line_length = len(i.rstrip())
            index += len(i)
            length += line_length

            if is_populated:
                continue

            is_populated = True
            result[current_chromosome]["line_length"] = line_length

        result[current_chromosome]["length"] = length

        for i in result:
            chromosome = Chromosome(
                self.fasta_file,
                name=i,
                line_length=result[i]["line_length"],
                label_length=result[i]["label_length"],
                index_start=result[i]["index_start"],
                length=result[i]["length"],
                labels=result[i]["labels"],
            )

            self.chromosomes[i] = chromosome

        return result

    def __getitem__(self, key):
        assert isinstance(key, str)

        return self.chromosomes[self._labels[key]]

    def sequence(
        self,
        chromosome: str,
        nucleotide: str,
        pos: int,
        from_nuc: int,
        to_nuc: int,
    ) -> list:
        """Gets a sequence from the fasta file by the position of a nucleotide on a
        chromosome, with a specified prefix length and suffix length.

        Parameters
        ----------
        chromosome : str
            The chromosome where the sequence is going to be obtained.
        nucleotide : str
            Label of the nucleotide.
        pos : int
            Position of the nucleotide on the chromosome.
        from_nuc : int
            Length of the prefix of the sequence.
        to_nuc : int
            Length of the suffix of the sequence.

        Returns
        -------
        The sequence divided in (prefix, nucleotide, suffix).
        """
        return self[chromosome].sequence(nucleotide, pos, from_nuc, to_nuc)
