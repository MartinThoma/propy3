# -*- coding: utf-8 -*-
"""
This module is used for computing the quasi sequence order descriptors based on
the given protein sequence. We can obtain two types of descriptors:
Sequence-order-coupling number and quasi-sequence-order descriptors. Two
distance matrixes between 20 amino acids are employed.

References
----------
.. [1] Kuo-Chen Chou. Prediction of Protein Subcellar Locations by Incorporating
       Quasi-Sequence-Order Effect. Biochemical and Biophysical Research
       Communications 2000, 278, 477-483.

.. [2] Kuo-Chen Chou and Yu-Dong Cai. Prediction of Protein sucellular locations by
       GO-FunD-PseAA predictor, Biochemical and Biophysical Research
       Communications, 2004, 320, 1236-1239.

.. [3] Gisbert Schneider and Paul wrede. The Rational Design of Amino Acid
       Sequences by Artifical Neural Networks and Simulated Molecular
       Evolution: Do Novo Design of an Idealized Leader Cleavge Site. Biophys
       Journal, 1994, 66, 335-344.

Authors: Dongsheng Cao and Yizeng Liang.
Date: 2012.09.03
Email: oriental-cds@163.com
"""

# Core Library
import json
import math
from typing import Any, Dict

# Third party
from pkg_resources import resource_filename

AALetter = list("ARNDCEQGHILKMFPSTWYV")

# Distance is the Schneider-Wrede physicochemical distance matrix
# used by Chou et. al.
filepath = resource_filename(
    __name__, "data/schneider-wrede-physicochemical-distance-matrix.json"
)
with open(filepath, "r") as f:
    _Distance1: Dict[str, float] = json.load(f)

# Distance is the Grantham chemical distance matrix used by Grantham et. al.
filepath = resource_filename(__name__, "data/grantham-chemical-distance-matrix.json")
with open(filepath, "r") as f:
    _Distance2: Dict[str, float] = json.load(f)


def GetSequenceOrderCouplingNumber(ProteinSequence, d=1, distancematrix=_Distance1):
    """
    Computing the dth-rank sequence order coupling number for a protein.

    Parameters
    ----------
    ProteinSequence : a pure protein sequence
    d : the gap between two amino acids.

    Returns
    -------
    a numeric value

    Example
    -------
    >>> result = GetSequenceOrderCouplingNumber(protein,d)
    """
    NumProtein = len(ProteinSequence)
    tau = 0.0
    for i in range(NumProtein - d):
        temp1 = ProteinSequence[i]
        temp2 = ProteinSequence[i + d]
        tau = tau + math.pow(distancematrix[temp1 + temp2], 2)
    return round(tau, 3)


def GetSequenceOrderCouplingNumberp(ProteinSequence, maxlag=30, distancematrix={}):
    """
    Computing the sequence order coupling numbers from 1 to maxlag

    for a given protein sequence based on the user-defined property.

    Parameters
    ----------
    protein : a pure protein sequence
    maxlag : int, optional (default: 30)
        the maximum lag and the length of the protein should be larger
        than maxlag.
    distancematrix : Dict
        a dict form containing 400 distance values

    Returns
    -------
    a dict form containing all sequence order coupling numbers based on the given property

    Examples
    --------
    >>> result = GetSequenceOrderCouplingNumberp(protein, maxlag,distancematrix)
    """
    # NumProtein = len(ProteinSequence) # TODO: this was calculated, but not
    # used... is here a bug?
    Tau = {}
    for i in range(maxlag):
        Tau["tau" + str(i + 1)] = GetSequenceOrderCouplingNumber(
            ProteinSequence, i + 1, distancematrix
        )
    return Tau


def GetSequenceOrderCouplingNumberSW(
    ProteinSequence, maxlag=30, distancematrix=_Distance1
):
    """
    Computing the sequence order coupling numbers from 1 to maxlag
    for a given protein sequence based on the Schneider-Wrede physicochemical
    distance matrix

    Examples
    --------
    >>> result = GetSequenceOrderCouplingNumberSW(protein, maxlag, distancematrix)

    Parameters
    ----------
    protein is a pure protein sequence

    maxlag is the maximum lag and the length of the protein should be larger
    than maxlag. default is 30.

    distancematrix is a dict form containing Schneider-Wrede physicochemical
    distance matrix. omitted!

    Returns
    -------
    result is a dict form containing all sequence order coupling
    numbers based on the Schneider-Wrede physicochemical distance matrix
    """
    # NumProtein = len(ProteinSequence)  # TODO: This was calculated, but not
    # used ... is here a bug?
    Tau = {}
    for i in range(maxlag):
        Tau["tausw" + str(i + 1)] = GetSequenceOrderCouplingNumber(
            ProteinSequence, i + 1, distancematrix
        )
    return Tau


def GetSequenceOrderCouplingNumberGrant(
    ProteinSequence, maxlag=30, distancematrix=_Distance2
):
    """
    Computing the sequence order coupling numbers from 1 to maxlag for a given
    protein sequence based on the Grantham chemical distance matrix.

    Examples
    --------
    >>> result = GetSequenceOrderCouplingNumberGrant(protein, maxlag, distancematrix)

    Parameters
    ----------
    protein is a pure protein sequence

    maxlag is the maximum lag and the length of the protein should be larger

    than maxlag. default is 30.

    distancematrix is a dict form containing Grantham chemical distance
    matrix. omitted!

    Returns
    -------
    result is a dict form containing all sequence order coupling numbers
    based on the Grantham chemical distance matrix
    """
    # NumProtein = len(ProteinSequence)  # TODO: This was calculated, but not
    # used ... is here a bug?
    Tau = {}
    for i in range(maxlag):
        Tau["taugrant" + str(i + 1)] = GetSequenceOrderCouplingNumber(
            ProteinSequence, i + 1, distancematrix
        )
    return Tau


def GetSequenceOrderCouplingNumberTotal(
    ProteinSequence, maxlag: int = 30
) -> Dict[Any, Any]:
    """
    Computing the sequence order coupling numbers from 1 to maxlag for a given
    protein sequence.

    Parameters
    ----------
    ProteinSequence : a pure protein sequence
    maxlag : int, optional (default: 30)
        the maximum lag and the length of the protein should be larger

    Returns
    -------
    result : Dict
        contains all sequence order coupling numbers

    Examples
    --------
    >>> result = GetSequenceOrderCouplingNumberTotal(protein, maxlag)
    """
    Tau: Dict[Any, Any] = {}
    Tau.update(GetSequenceOrderCouplingNumberSW(ProteinSequence, maxlag=maxlag))
    Tau.update(GetSequenceOrderCouplingNumberGrant(ProteinSequence, maxlag=maxlag))
    return Tau


def GetAAComposition(ProteinSequence) -> Dict[str, float]:
    """
    Calculate the composition of Amino acids for a given protein sequence.

    Parameters
    ----------
    ProteinSequence : a pure protein sequence.

    Returns
    -------
    result : Dict[str, float]
        contains the composition of 20 amino acids.

    Examples
    --------
    >>> result = CalculateAAComposition(protein)
    """
    LengthSequence = len(ProteinSequence)
    Result: Dict[str, float] = {}
    for i in AALetter:
        Result[i] = round(float(ProteinSequence.count(i)) / LengthSequence, 3)
    return Result


def GetQuasiSequenceOrder1(
    ProteinSequence, maxlag: int = 30, weight: float = 0.1, distancematrix={}
):
    """
    Computing the first 20 quasi-sequence-order descriptors for a given protein
    sequence.

    Examples
    --------
    >>> result = GetQuasiSequenceOrder1(protein, maxlag, weigt)

    see method GetQuasiSequenceOrder for the choice of parameters.
    """
    rightpart = 0.0
    for i in range(maxlag):
        rightpart = rightpart + GetSequenceOrderCouplingNumber(
            ProteinSequence, i + 1, distancematrix
        )
    AAC = GetAAComposition(ProteinSequence)
    result: Dict[str, float] = {}
    temp = 1 + weight * rightpart
    for index, aaletter_char in enumerate(AALetter):
        result["QSO" + str(index + 1)] = round(AAC[aaletter_char] / temp, 6)

    return result


def GetQuasiSequenceOrder2(ProteinSequence, maxlag=30, weight=0.1, distancematrix={}):
    """
    Computing the last maxlag quasi-sequence-order descriptors for a given
    protein sequence.

    Examples
    --------
    >>> result = GetQuasiSequenceOrder2(protein,maxlag,weigt)

    see method GetQuasiSequenceOrder for the choice of parameters.
    """
    rightpart = []
    for i in range(maxlag):
        rightpart.append(
            GetSequenceOrderCouplingNumber(ProteinSequence, i + 1, distancematrix)
        )
    # AAC = GetAAComposition(ProteinSequence)  # TODO: was not used. Bug?
    result = {}
    temp = 1 + weight * sum(rightpart)
    for index in range(20, 20 + maxlag):
        result["QSO" + str(index + 1)] = round(weight * rightpart[index - 20] / temp, 6)
    return result


def GetQuasiSequenceOrder1SW(
    ProteinSequence, maxlag=30, weight=0.1, distancematrix=_Distance1
):
    """
    Computing the first 20 quasi-sequence-order descriptors for a given protein
    sequence.

    Examples
    --------
    >>> result = GetQuasiSequenceOrder1SW(protein, maxlag, weigt)

    see method GetQuasiSequenceOrder for the choice of parameters.
    """
    rightpart = 0.0
    for i in range(maxlag):
        rightpart = rightpart + GetSequenceOrderCouplingNumber(
            ProteinSequence, i + 1, distancematrix
        )
    AAC = GetAAComposition(ProteinSequence)
    result = {}
    temp = 1 + weight * rightpart
    for index, aaletter_char in enumerate(AALetter):
        result["QSOSW" + str(index + 1)] = round(AAC[aaletter_char] / temp, 6)

    return result


def GetQuasiSequenceOrder2SW(
    ProteinSequence, maxlag=30, weight=0.1, distancematrix=_Distance1
):
    """
    Computing the last maxlag quasi-sequence-order descriptors for
    a given protein sequence.

    Examples
    --------
    >>> result = GetQuasiSequenceOrder2SW(protein,maxlag,weigt)

    see method GetQuasiSequenceOrder for the choice of parameters.
    """
    rightpart = []
    for i in range(maxlag):
        rightpart.append(
            GetSequenceOrderCouplingNumber(ProteinSequence, i + 1, distancematrix)
        )
    # AAC = GetAAComposition(ProteinSequence)  # TODO: was not used. Bug?
    result = {}
    temp = 1 + weight * sum(rightpart)
    for index in range(20, 20 + maxlag):
        result["QSOSW" + str(index + 1)] = round(
            weight * rightpart[index - 20] / temp, 6
        )

    return result


def GetQuasiSequenceOrder1Grant(
    ProteinSequence, maxlag=30, weight=0.1, distancematrix=_Distance2
):
    """
    Computing the first 20 quasi-sequence-order descriptors for
    a given protein sequence.

    Examples
    --------
    >>> result = GetQuasiSequenceOrder1Grant(protein, maxlag, weigt)

    see method GetQuasiSequenceOrder for the choice of parameters.
    """
    rightpart = 0.0
    for i in range(maxlag):
        rightpart = rightpart + GetSequenceOrderCouplingNumber(
            ProteinSequence, i + 1, distancematrix
        )
    AAC = GetAAComposition(ProteinSequence)
    result = {}
    temp = 1 + weight * rightpart
    for index, aaletter_char in enumerate(AALetter):
        result["QSOgrant" + str(index + 1)] = round(AAC[aaletter_char] / temp, 6)

    return result


def GetQuasiSequenceOrder2Grant(
    ProteinSequence, maxlag=30, weight=0.1, distancematrix=_Distance2
):
    """
    Computing the last maxlag quasi-sequence-order descriptors for
    a given protein sequence.

    Examples
    --------
    >>> result = GetQuasiSequenceOrder2Grant(protein, maxlag, weigt)

    see method GetQuasiSequenceOrder for the choice of parameters.
    """
    rightpart = []
    for i in range(maxlag):
        rightpart.append(
            GetSequenceOrderCouplingNumber(ProteinSequence, i + 1, distancematrix)
        )
    # AAC = GetAAComposition(ProteinSequence)  # TODO: Was not used. Bug?
    result = {}
    temp = 1 + weight * sum(rightpart)
    for index in range(20, 20 + maxlag):
        result["QSOgrant" + str(index + 1)] = round(
            weight * rightpart[index - 20] / temp, 6
        )

    return result


def GetQuasiSequenceOrder(ProteinSequence, maxlag=30, weight=0.1) -> Dict[Any, Any]:
    """
    Computing quasi-sequence-order descriptors for a given protein.

    Parameters
    ----------
    protein is a pure protein sequence

    maxlag : int, optional (default: 30)
        the maximum lag and the length of the protein should be larger than
        maxlag
    weight : float, optional (default: 0.1)
        a weight factor. Please see reference 1 for its choice.

    Returns
    -------
    result : Dict
        contains all quasi-sequence-order descriptors

    Examples
    --------
    >>> result = GetQuasiSequenceOrder(protein, maxlag, weight)

    References
    ----------
    .. [1] Kuo-Chen Chou. Prediction of Protein Subcellar Locations by
           Incorporating Quasi-Sequence-Order Effect. Biochemical and
           Biophysical Research Communications 2000, 278, 477-483.
    """
    result: Dict[Any, Any] = dict()
    result.update(GetQuasiSequenceOrder1SW(ProteinSequence, maxlag, weight, _Distance1))
    result.update(GetQuasiSequenceOrder2SW(ProteinSequence, maxlag, weight, _Distance1))
    result.update(
        GetQuasiSequenceOrder1Grant(ProteinSequence, maxlag, weight, _Distance2)
    )
    result.update(
        GetQuasiSequenceOrder2Grant(ProteinSequence, maxlag, weight, _Distance2)
    )
    return result


def GetQuasiSequenceOrderp(
    ProteinSequence,
    maxlag: int = 30,
    weight: float = 0.1,
    distancematrix: Dict[Any, Any] = {},
):
    """
    Computing quasi-sequence-order descriptors for a given protein.

    Parameters
    ----------
    ProteinSequence : a pure protein sequence
    maxlag : int, optional (default: 30)
        the maximum lag and the length of the protein should be larger than
        maxlag
    weight : float, optional (default: 0.1)
        a weight factor. Please see reference 1 for its choice.
    distancematrix : Dict[Any, Any]
        contains 400 distance values

    Returns
    -------
    result is a dict form containing all quasi-sequence-order descriptors

    Examples
    --------
    >>> result = GetQuasiSequenceOrderp(protein,maxlag,weight,distancematrix)

    References
    ----------
    .. [1] Kuo-Chen Chou. Prediction of Protein Subcellar Locations by
           Incorporating Quasi-Sequence-Order Effect. Biochemical and
           Biophysical Research Communications 2000, 278, 477-483.
    """
    result: Dict[Any, Any] = dict()
    result.update(
        GetQuasiSequenceOrder1(ProteinSequence, maxlag, weight, distancematrix)
    )
    result.update(
        GetQuasiSequenceOrder2(ProteinSequence, maxlag, weight, distancematrix)
    )
    return result
