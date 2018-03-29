def read_sequences(path):
    """
    Read common text file from given path and return generator with all sequences from file
    :param path: str - path to file
    :return: generator - generator with sequences
    """
    with open(path) as source:
        return (x.strip() for x in source.readlines())