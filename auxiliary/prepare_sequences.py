from Bio import SeqIO


def prepare_sequences(path_to_first_file, path_to_second_file):
    """
    Read fasta files and return them. This function is made specially for align_one_vs_all
    :param path_to_first_file: str - path to fasta file, by convention let`s supply here templates
    :param path_to_second_file: str - path to fasta file, by convention let`s supply here reads
    :return: tuple - generators with SeqRecord objects
    """
    # Read files
    first = SeqIO.parse(path_to_first_file, 'fasta')
    second = SeqIO.parse(path_to_second_file, 'fasta')
    return first, second