# -*- coding: utf-8 -*-
"""
The prediction of functional sites (e.g.,methylation) of proteins usually needs
to split the total protein into a set of segments around specific amino acid.
Given a specific window size p, we can obtain all segments of length equal to
(2*p+1) very easily. Note that the output of the method is a list form.

Authors: Dongsheng Cao and Yizeng Liang.
Date: 2012.9.2
Email: oriental-cds@163.com
"""

# Core Library
import re
from typing import List

AALetter = [
    "A",
    "R",
    "N",
    "D",
    "C",
    "E",
    "Q",
    "G",
    "H",
    "I",
    "L",
    "K",
    "M",
    "F",
    "P",
    "S",
    "T",
    "W",
    "Y",
    "V",
]


def GetSubSequence(ProteinSequence: str, ToAA: str = "S", window: int = 3) -> List[str]:
    """
    Get all 2*window+1 sub-sequences whose cener is ToAA in a protein.

    Parameters
    ----------
    ProteinSequence : str
        a pure problem sequence
    ToAA :str
        the central (query point) amino acid in the sub-sequence
    window : int
        the span

    Returns
    -------
    result : List[str]
        contains all satisfied sub-sequences

    Examples
    --------
    >>> result = GetSubSequence(protein, ToAA, window)
    """

    if ToAA not in AALetter:
        ToAA = ProteinSequence[1]

    Num = len(ProteinSequence)
    seqiter = re.finditer(ToAA, ProteinSequence)
    AAindex: List[int] = []
    for seq_element in seqiter:
        AAindex.append(seq_element.end())

    result = []
    for i in AAindex:
        if i - window > 0 and Num - i + 1 - window > 0:
            temp = ProteinSequence[i - window - 1 : i + window]
            result.append(temp)

    return result
