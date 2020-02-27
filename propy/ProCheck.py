# -*- coding: utf-8 -*-
"""
Check whether the input protein sequence is a valid amino acid sequence.

Authors: Dongsheng Cao and Yizeng Liang.
Date: 2012.09.09
Email: oriental-cds@163.com
"""

AALetter = list("ARNDCEQGHILKMFPSTWYV")


def ProteinCheck(ProteinSequence: str) -> int:
    """
    Check whether the protein sequence is a valid amino acid sequence or not.

    Parameters
    ----------
    ProteinSequence : a pure protein sequence

    Returns
    -------
    flag : bool
        if the check is no problem, result will return the length of protein.
        if the check has problems, result will return 0.

    Examples
    --------
    >>> result = ProteinCheck(protein)
    """
    NumPro = len(ProteinSequence)
    for i in ProteinSequence:
        if i not in AALetter:
            flag = 0
            break
        else:
            flag = NumPro

    return flag
