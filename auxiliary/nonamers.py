from auxiliary.read_sequences import read_sequences


def nonamers(path_to_nonamers='../data/conserve/hv9'):
    """
    Read common text file from given path and return generator with all sequences from file
    :param path_to_nonamers: str - path to file
    :return: generator - generator with sequences
    """
    return read_sequences(path_to_nonamers)