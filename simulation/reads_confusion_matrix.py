from Bio import pairwise2
from simulation.categorize_reads import categorize_reads


# Use categorize reads 2 times and get confusion matrix for all threshold levels

def confusion_matrix(overlapped='/home/arleg/ig_construction/alignment/overlapped',
                     nonoverlapped='/home/arleg/ig_construction/alignment/non_overlapped',
                     genes='/home/arleg/ig_construction/alignment/test_genes.fa',
                     ceil=101,
                     score_function=pairwise2.align.localms,
                     *args, **kwargs):
    """
    Get basic parameters - numbers of true positive, false negative, false positive and true negative
    :param overlapped: str - path to file with overlapped reads
    :param nonoverlapped: str - path to file with nonoverlapped reads
    :param genes: str - path to file with genes
    :param ceil: int - upper boundary of threshold in testing in truthness of classification
    :param score_function: function - function which will make and score alignment
    :param args: unnamed args to function
    :param kwargs: named args to function
    :return: tuple - 4 np arrays (, ceil) with tp, fn, fp, tn
    """
    # Find 4 core characteristics via categorize_reads function upon overlapped and non-overlapped reads
    tp, fn = categorize_reads(overlapped, genes, True, ceil, score_function, *args, **kwargs)
    tn, fp = categorize_reads(nonoverlapped, genes, False, ceil, score_function, *args, **kwargs)

    return tp, fn, fp, tn


# a = confusion_matrix('/home/arleg/ig_construction/alignment/test_overlap',
#                  '/home/arleg/ig_construction/alignment/test_non_overlap',
#                  '/home/arleg/ig_construction/alignment/test_genes.fa',
#                  101,
#                  pairwise2.align.localms, 1, -2, -5, -1, score_only=True)
# for i in a:
#     print("Characteristic:", i)