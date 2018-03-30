def precision(tp, fp):
    """
    Compute precision given true positives and false positives
    :param tp: np array - true positives
    :param fp: np array - false positives
    :return: np array - precisions
    """
    return tp / (tp + fp)
