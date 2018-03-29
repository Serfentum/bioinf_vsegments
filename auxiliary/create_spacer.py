import random


def create_spacer(length, monomers=('A', 'T', 'G', 'C'), cum_distribution=(0.25, 0.5, 0.75, 1)):
    """
    Create random polymer with specified length from given monomers with given distribution
    :param length: int - polymer length
    :param monomers: sequence - set of monomers to built polymer from
    :param cum_distribution: sequence - according cumulative probabilities for monomers
    :return: str - polymer sequence
    """
    return ''.join(random.choices(monomers, cum_weights=cum_distribution, k=length))