from Bio import SeqIO


def get_sequences(path_to_genes='../data/main/simple_fasta/all/ig_hv_all'):
    """
    Read fasta from given path and return generator with all sequences from file
    :param path_to_genes: str - path to fasta file
    :return: generator - generator with sequences
    """
    return (str(x.seq) for x in SeqIO.parse(path_to_genes, 'fasta'))