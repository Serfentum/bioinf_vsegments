import random


def sequencing(sequence, read_length, error_rate=0.02):
    """
    Return all reads from given sequence with some errors during sequencing
    :param sequence: str - sequence of polymer
    :param read_length: int - length of reads
    :param error_rate: float - fraction of errors during sequencing
    :return: list - list of reads
    """
    # Generate slices of original sequence with length equal to read_length where each nucleotide can be mutated
    return [sequence_read(sequence, start, read_length, error_rate)
            for start in range(len(sequence) - read_length + 1)]


def sequence_read(sequence, start, read_length, error_rate=0.02):
    """
    Return read sequenced from sequence with some errors during sequencing
    :param sequence: str - sequence of polymer
    :param start: int - index from which position sequencing
    :param read_length: int - length of read
    :param error_rate: float - fraction of errors during sequencing
    :return: str - read
    """
    return ''.join([nucl
                    if random.random() > error_rate
                    else random.choice(list(filter(lambda x: x != nucl, ['A', 'T', 'G', 'C'])))
                    for nucl in sequence[start:start + read_length]])