import random
from Bio import SeqIO


def create_spacer(length, monomers=('A', 'T', 'G', 'C'), cum_distribution=(0.25, 0.5, 0.75, 1)):
    """
    Create random polymer with specified length from given monomers with given distribution
    :param length: int - polymer length
    :param monomers: sequence - set of monomers to built polymer from
    :param cum_distribution: sequence - according cumulative probabilities for monomers
    :return: str - polymer sequence
    """
    return ''.join(random.choices(monomers, cum_weights=cum_distribution, k=length))


def get_sequences(path_to_genes='../data/main/simple_fasta/all/ig_hv_all'):
    """
    Read fasta from given path and return generator with all sequences from file
    :param path_to_genes: str - path to fasta file
    :return: generator - generator with sequences
    """
    return (x.seq for x in SeqIO.parse(path_to_genes, 'fasta'))


def read_sequences(path):
    """
    Read common text file from given path and return generator with all sequences from file
    :param path: str - path to file
    :return: generator - generator with sequences
    """
    with open(path) as source:
        return (x.strip() for x in source.readlines())


def heptamers(path_to_heptamers='../data/conserve/hv7'):
    """
    Read common text file from given path and return generator with all sequences from file
    :param path_to_heptamers: str - path to file
    :return: generator - generator with sequences
    """
    return read_sequences(path_to_heptamers)


def nonamers(path_to_nonamers='../data/conserve/hv9'):
    """
    Read common text file from given path and return generator with all sequences from file
    :param path_to_nonamers: str - path to file
    :return: generator - generator with sequences
    """
    return read_sequences(path_to_nonamers)
