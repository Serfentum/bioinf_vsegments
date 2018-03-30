import numpy as np
import matplotlib.pyplot as plt
from simulation.recall import recall
from simulation.false_positive_rate import false_positive_rate


# For testing
# tp = np.fromstring('658 658 658 658 658 658 657 574 368 223 102  52  29  21  17  12   9   5  5   3   3   3   3   3   2   2   2   2   2   2   2   2   2   2   2   2  2   2   2   2   1   1   1   1   1   1   1   1   1   1   1   1   1   1 1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   0   0  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  0   0   0   0   0   0   0   0   0   0   0', sep=' ')
# fn = np.fromstring('''0   0   0   0   0   0   1  84 290 435 556 606 629 637 641 646 649 653
#  653 655 655 655 655 655 656 656 656 656 656 656 656 656 656 656 656 656
#  656 656 656 656 657 657 657 657 657 657 657 657 657 657 657 657 657 657
#  657 657 657 657 657 657 657 657 657 657 657 657 657 657 657 657 658 658
#  658 658 658 658 658 658 658 658 658 658 658 658 658 658 658 658 658 658
#  658 658 658 658 658 658 658 658 658 658 658''', sep=' ')
# fp = np.fromstring('''64 64 64 64 64 64 60 43 25  7  4  0  0  0  0  0  0  0  0  0  0  0  0  0
#   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#   0  0  0  0  0''', sep=' ')
# tn = np.fromstring(''' 0  0  0  0  0  0  4 21 39 57 60 64 64 64 64 64 64 64 64 64 64 64 64 64
#  64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64
#  64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64
#  64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64
#  64 64 64 64 64''', sep=' ')


def plot_roc(tp, fp, fn, tn, output):
    """
    Draw ROC curve and no-discrimination line with grid from main parameters
    :param tp: np array - true positives
    :param fp: np array - false positives
    :param fn: np array - false negatives
    :param tn: np array - true negatives
    :param output: str - path to output file with its name
    :return:
    """
    # Computing false positive rates and true positive rates
    fpr = false_positive_rate(tn, fp)
    tpr = recall(tp, fn)

    # Plotting
    plt.plot(fpr, tpr)
    plt.plot([0, 1], [0, 1])
    plt.grid()

    # Decorate plot
    plt.title('ROC')
    plt.xlabel('False positive rate')
    plt.ylabel('True positive rate')

    # Save plot
    plt.savefig(output)


# plot_roc(tp, fp, fn, tn, 'pr.svg')

