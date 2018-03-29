import numpy as np
from Bio import pairwise2
from auxiliary.get_sequences import get_sequences
from auxiliary.read_sequences import read_sequences


# This function is made to get data for confusion matrix. Reads are scored by scoring function, after this they are
# compared with threshold and counting. Threshold ranges from 0 to ceil inclusively.
# Each threshold will be its own function on ROC and numbers of passed and unpassed reads - tpr and fpr rates
# correspondingly.
# It is possible to specify another scoring function except simple alignment score like here, yet it should make an
# alignment first.


def categorize_reads(reads='/home/arleg/ig_construction/alignment/overlapped',
                     genes='/home/arleg/ig_construction/alignment/test_genes.fa',
                     examples_true=True,
                     ceil=101,
                     score_function=pairwise2.align.localms,
                     *args,
                     **kwargs):

    """
    Count reads from file
    :param reads: str - path to file with reads
    :param genes: str - path to file with genes
    :param examples_true: boolean - whether reads was from region or not
    :param ceil: int - upper boundary of threshold in testing in truthness of classification
    :param score_function: function - function which will make and score alignment
    :param args: unnamed args to function
    :param kwargs: named args to function
    :return: tuple - 2 np arrays (, ceil) with number of reads scored greater than threshold or lesser on each threshold
    from 0 to ceil including; 1st going array with numbers of true answers
    """
    # Initialize containers for output
    mapped = np.zeros(ceil, int)
    unmapped = np.zeros(ceil, int)

    # Read sequences of genes and reads
    genes = get_sequences(genes)
    reads = read_sequences(reads)

    # Align all genes vs all reads and score alignment.
    # Than allocate this read to mapped or unmapped for all thresholds depending on comparison of score and threshold
    for gene in genes:
        for read in reads:
            score = score_function(read, gene, *args, **kwargs)
            for threshold in range(ceil):
                if score >= threshold:
                    mapped[threshold] += 1
                else:
                    unmapped[threshold] += 1

    # If passed to function reads was really mapped than mapped is a true positives and unmapped is a false negatives,
    # otherwise unmapped is a true negatives and mapped is a false positives
    if examples_true:
        return mapped, unmapped
    return unmapped, mapped


# a = categorize_reads('/home/arleg/ig_construction/alignment/test_overlap',
#                      '/home/arleg/ig_construction/alignment/test_genes.fa',
#                      True,
#                      101,
#                      pairwise2.align.localms,
#                      1, -2, -5, -1, score_only=True)

# for i in a:
#     print(i)




