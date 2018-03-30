def recall(tp, fn):
    """
    Compute recall given true positives and false negatives.
    Recall is the same characteristic as true positive rate
    :param tp: np array - true positives
    :param fn: np array - false negatives
    :return: np array - recalls
    """
    return tp / (tp + fn)
