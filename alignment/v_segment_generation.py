from auxiliary import *


def combinations(path_to_genes, path_to_heptamers, path_to_nonamers, length, monomers, cum_distribution):
    """
    Make all combinations of given genes, heptamers, random spacer and nonamers
    Template - {core}{heptamer}{random_spacer}{nonamer}
    :param path_to_genes: str - path to fasta file
    :param path_to_heptamers: str - path to file
    :param path_to_nonamers: str - path to file
    :param length: int - polymer length
    :param monomers: sequence - set of monomers to built polymer from
    :param cum_distribution: sequence - according cumulative probabilities for monomers
    :return:
    """
    # Combine sequences from gene_cores, heptamers, create_spacer and nonamer into V-segment
    combs = (f'{core}{heptamer}{create_spacer(length, monomers, cum_distribution)}{nonamer}'
                    for core in gene_cores(path_to_genes)
                    for heptamer in heptamers(path_to_heptamers)
                    for nonamer in nonamers(path_to_nonamers))
    return combs


def v_segments():
    """
    Default version of combinations() to make full V segments
    """
    return combinations(path_to_genes='../data/main/simple_fasta/all/ig_hv_all',
                        path_to_heptamers='../data/conserve/hv7',
                        path_to_nonamers='../data/conserve/hv9',
                        length=23,
                        monomers=('A', 'T', 'G', 'C'),
                        cum_distribution=(0.25, 0.5, 0.75, 1))

# path_to_heptamers='../data/conserve/hv7'
# path_to_nonamers='../data/conserve/hv9'
# monomers=('A', 'T', 'G', 'C')
# cum_distribution=(0.25, 0.5, 0.75, 1)
# for i, s in enumerate(v_segments()):
#     print(i, s)
#     if i > 30:
#         break
