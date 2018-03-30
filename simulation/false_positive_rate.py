def false_positive_rate(tn, fp):
    """
    Compute false positive rate given true negatives and false positives
    :param tn: np array - true negatives
    :param fp: np array - alse positives
    :return: np array - false positive rates
    """
    return fp / (tn + fp)
