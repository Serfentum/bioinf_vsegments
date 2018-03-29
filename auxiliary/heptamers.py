from auxiliary.read_sequences import read_sequences


def heptamers(path_to_heptamers='../data/conserve/hv7'):
    """
    Read common text file from given path and return generator with all sequences from file
    :param path_to_heptamers: str - path to file
    :return: generator - generator with sequences
    """
    return read_sequences(path_to_heptamers)