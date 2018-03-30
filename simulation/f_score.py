def f_score(precision, recall, beta=1):
    """
    Compute F-score with specified parameter
    :param precision: np array - 1d array with precision
    :param recall: np array - 1d array with recall
    :param beta: float - value of beta parameter - weight of precision in scoring
    :return: float - F-score
    """
    factor = pow(beta, 2)
    return (1 + factor) * precision * recall / (factor * precision + recall)


# Testing
# from simulation.precision import precision
# from simulation.recall import recall
# import numpy as np
#
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
#
# p = precision(tp, fp)
# r = recall(tp, fn)
#
# print(f_score(p, r, 0.5))
# print(f_score(p, r))
# print(f_score(p, r, 2))