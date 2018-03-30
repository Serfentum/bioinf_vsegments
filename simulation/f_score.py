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