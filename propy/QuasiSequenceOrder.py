# propy3, formerly protpy, is a Python package to compute protein descriptors
# Copyright (C) 2012 Dongsheng Cao and Yizeng Liang, oriental-cds@163.com
# Copyright (C) 2020-2022 Martin Thoma

# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; in version 2
# of the License.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor,
# Boston, MA  02110-1301, USA.
"""
Compute the quasi sequence order descriptors based on the given protein
sequence. We can obtain two types of descriptors: Sequence-order-coupling
number and quasi-sequence-order descriptors. Two distance matrixes between 20
amino acids are employed.

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
"""

# Core Library
import json
import math
from typing import Any, Dict

# Third party
from pkg_resources import resource_filename

# First party
from propy import AALetter

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
    _Distance2: Dict[str, int] = json.load(f)


def GetSequenceOrderCouplingNumber(
    ProteinSequence: str, d: int = 1, distancematrix: Dict[str, float] = _Distance1
):
    """
    Compute the dth-rank sequence order coupling number for a protein.

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence
    d : int
        the gap between two amino acids.
    distancematrix : Dict[str, float]

    Returns
    -------
    tau : float

    Example
    -------
    >>> from propy.GetProteinFromUniprot import GetProteinSequence
    >>> protein = GetProteinSequence(ProteinID="Q9NQ39")
    >>> result = GetSequenceOrderCouplingNumber(protein)
    """
    NumProtein = len(ProteinSequence)
    tau = 0.0
    for i in range(NumProtein - d):
        temp1 = ProteinSequence[i]
        temp2 = ProteinSequence[i + d]
        tau = tau + math.pow(distancematrix[temp1 + temp2], 2)
    return round(tau, 3)


def GetSequenceOrderCouplingNumberp(
    ProteinSequence: str, maxlag: int = 30, distancematrix: Dict[Any, Any] = None
):
    """
    Compute the sequence order coupling numbers from 1 to maxlag
    for a given protein sequence based on the user-defined property.

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence
    maxlag : int, optional (default: 30)
        the maximum lag and the length of the protein should be larger
        than maxlag.
    distancematrix : Dict[Any, Any]
        contains 400 distance values

    Returns
    -------
    Tau : Dict[str]
        contains all sequence order coupling numbers based on the given property

    Examples
    --------
    >>> from propy.GetProteinFromUniprot import GetProteinSequence
    >>> protein = GetProteinSequence(ProteinID="Q9NQ39")
    >>> result = GetSequenceOrderCouplingNumberp(protein)
    """
    if distancematrix is None:
        distancematrix = {}
    Tau = {}
    for i in range(maxlag):
        Tau["tau" + str(i + 1)] = GetSequenceOrderCouplingNumber(
            ProteinSequence, i + 1, distancematrix
        )
    return Tau


def GetSequenceOrderCouplingNumberSW(
    ProteinSequence: str, maxlag: int = 30, distancematrix=_Distance1
):
    """
    Compute the sequence order coupling numbers from 1 to maxlag for a given
    protein sequence based on the Schneider-Wrede physicochemical distance
    matrix.

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence
    maxlag : int, optional (default: 30)
        the maximum lag and the length of the protein should be larger than
        maxlag
    distancematrix : Dict[Any, Any]
        contains Schneider-Wrede physicochemical distance matrix

    Returns
    -------
    Tau : Dict[Any, Any]
        contains all sequence order coupling numbers based on the
        Schneider-Wrede physicochemical distance matrix

    Examples
    --------
    >>> from propy.GetProteinFromUniprot import GetProteinSequence
    >>> protein = GetProteinSequence(ProteinID="Q9NQ39")
    >>> result = GetSequenceOrderCouplingNumberSW(protein)
    """
    Tau = {}
    for i in range(maxlag):
        Tau["tausw" + str(i + 1)] = GetSequenceOrderCouplingNumber(
            ProteinSequence, i + 1, distancematrix
        )
    return Tau


def GetSequenceOrderCouplingNumberGrant(
    ProteinSequence: str, maxlag: int = 30, distancematrix=_Distance2
):
    """
    Compute the sequence order coupling numbers from 1 to maxlag for a given
    protein sequence based on the Grantham chemical distance matrix.

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence
    maxlag : int, optional (default: 30)
        the maximum lag and the length of the protein should be larger than
        maxlag
    distancematrix : Dict[Any, Any]
        contains Schneider-Wrede physicochemical distance matrix

    Returns
    -------
    Tau : Dict[Any, Any]
        contains all sequence order coupling numbers based on the Grantham
        chemical distance matrix

    Examples
    --------
    >>> from propy.GetProteinFromUniprot import GetProteinSequence
    >>> protein = GetProteinSequence(ProteinID="Q9NQ39")
    >>> result = GetSequenceOrderCouplingNumberGrant(protein)
    """
    Tau = {}
    for i in range(maxlag):
        Tau["taugrant" + str(i + 1)] = GetSequenceOrderCouplingNumber(
            ProteinSequence, i + 1, distancematrix
        )
    return Tau


def GetSequenceOrderCouplingNumberTotal(
    ProteinSequence: str, maxlag: int = 30
) -> Dict[Any, Any]:
    """
    Compute the sequence order coupling numbers from 1 to maxlag for a given
    protein sequence.

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence
    maxlag : int, optional (default: 30)
        the maximum lag and the length of the protein should be larger

    Returns
    -------
    result : Dict
        contains all sequence order coupling numbers

    Examples
    --------
    >>> from propy.GetProteinFromUniprot import GetProteinSequence
    >>> protein = GetProteinSequence(ProteinID="Q9NQ39")
    >>> result = GetSequenceOrderCouplingNumberTotal(protein)
    """
    Tau: Dict[Any, Any] = {}
    Tau.update(GetSequenceOrderCouplingNumberSW(ProteinSequence, maxlag=maxlag))
    Tau.update(GetSequenceOrderCouplingNumberGrant(ProteinSequence, maxlag=maxlag))
    return Tau


def GetAAComposition(ProteinSequence: str) -> Dict[str, float]:
    """
    Calculate the composition of Amino acids for a given protein sequence.

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    Returns
    -------
    result : Dict[str, float]
        contains the composition of 20 amino acids.

    Examples
    --------
    >>> from propy.GetProteinFromUniprot import GetProteinSequence
    >>> from propy.AAComposition import CalculateAAComposition
    >>> protein = GetProteinSequence(ProteinID="Q9NQ39")
    >>> result = CalculateAAComposition(protein)
    """
    LengthSequence = len(ProteinSequence)
    result: Dict[str, float] = {}
    for i in AALetter:
        result[i] = round(float(ProteinSequence.count(i)) / LengthSequence, 3)
    return result


def GetQuasiSequenceOrder1(
    ProteinSequence: str, maxlag: int = 30, weight: float = 0.1, distancematrix=None
):
    """
    Compute the first 20 quasi-sequence-order descriptors for a given protein
    sequence.

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence
    maxlag : int, optional (default: 30)
        the maximum lag and the length of the protein should be larger

    Examples
    --------
    >>> from propy.GetProteinFromUniprot import GetProteinSequence
    >>> protein = GetProteinSequence(ProteinID="Q9NQ39")
    >>> result = GetQuasiSequenceOrder1(protein)

    see :py:func:`GetQuasiSequenceOrder` for the choice of parameters.
    """
    if distancematrix is None:
        distancematrix = {}
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


def GetQuasiSequenceOrder2(
    ProteinSequence: str, maxlag=30, weight=0.1, distancematrix=None
):
    """
    Compute the last maxlag quasi-sequence-order descriptors for a given
    protein sequence.

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence
    maxlag : int, optional (default: 30)
        the maximum lag and the length of the protein should be larger

    Examples
    --------
    >>> from propy.GetProteinFromUniprot import GetProteinSequence
    >>> protein = GetProteinSequence(ProteinID="Q9NQ39")
    >>> result = GetQuasiSequenceOrder2(protein)

    see :py:func:`GetQuasiSequenceOrder` for the choice of parameters.
    """
    if distancematrix is None:
        distancematrix = {}
    rightpart = []
    for i in range(maxlag):
        rightpart.append(
            GetSequenceOrderCouplingNumber(ProteinSequence, i + 1, distancematrix)
        )
    result = {}
    temp = 1 + weight * sum(rightpart)
    for index in range(20, 20 + maxlag):
        result["QSO" + str(index + 1)] = round(weight * rightpart[index - 20] / temp, 6)
    return result


def GetQuasiSequenceOrder1SW(
    ProteinSequence: str, maxlag=30, weight=0.1, distancematrix=_Distance1
):
    """
    Compute the first 20 quasi-sequence-order descriptors for a given protein
    sequence.

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence
    maxlag : int, optional (default: 30)
        the maximum lag and the length of the protein should be larger

    Examples
    --------
    >>> from propy.GetProteinFromUniprot import GetProteinSequence
    >>> protein = GetProteinSequence(ProteinID="Q9NQ39")
    >>> result = GetQuasiSequenceOrder1SW(protein)

    see :py:func:`GetQuasiSequenceOrder` for the choice of parameters.
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
    ProteinSequence: str, maxlag=30, weight=0.1, distancematrix=_Distance1
):
    """
    Compute the last maxlag quasi-sequence-order descriptors for a given
    protein sequence.

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence
    maxlag : int, optional (default: 30)
        the maximum lag and the length of the protein should be larger

    Examples
    --------
    >>> from propy.GetProteinFromUniprot import GetProteinSequence
    >>> protein = GetProteinSequence(ProteinID="Q9NQ39")
    >>> result = GetQuasiSequenceOrder2SW(protein)

    see :py:func:`GetQuasiSequenceOrder` for the choice of parameters.
    """
    rightpart = []
    for i in range(maxlag):
        rightpart.append(
            GetSequenceOrderCouplingNumber(ProteinSequence, i + 1, distancematrix)
        )
    result = {}
    temp = 1 + weight * sum(rightpart)
    for index in range(20, 20 + maxlag):
        result["QSOSW" + str(index + 1)] = round(
            weight * rightpart[index - 20] / temp, 6
        )

    return result


def GetQuasiSequenceOrder1Grant(
    ProteinSequence: str,
    maxlag: int = 30,
    weight: float = 0.1,
    distancematrix=_Distance2,
):
    """
    Compute the first 20 quasi-sequence-order descriptors for a given protein
    sequence.

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence
    maxlag : int, optional (default: 30)
        the maximum lag and the length of the protein should be larger

    Examples
    --------
    >>> from propy.GetProteinFromUniprot import GetProteinSequence
    >>> protein = GetProteinSequence(ProteinID="Q9NQ39")
    >>> result = GetQuasiSequenceOrder1Grant(protein)

    see :py:func:`GetQuasiSequenceOrder` for the choice of parameters.
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
    ProteinSequence: str,
    maxlag: int = 30,
    weight: float = 0.1,
    distancematrix=_Distance2,
):
    """
    Compute the last maxlag quasi-sequence-order descriptors for a given
    protein sequence.

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence
    maxlag : int, optional (default: 30)
        the maximum lag and the length of the protein should be larger

    Examples
    --------
    >>> from propy.GetProteinFromUniprot import GetProteinSequence
    >>> protein = GetProteinSequence(ProteinID="Q9NQ39")
    >>> result = GetQuasiSequenceOrder2Grant(protein)

    see :py:func:`GetQuasiSequenceOrder` for the choice of parameters.
    """
    rightpart = []
    for i in range(maxlag):
        rightpart.append(
            GetSequenceOrderCouplingNumber(ProteinSequence, i + 1, distancematrix)
        )
    result = {}
    temp = 1 + weight * sum(rightpart)
    for index in range(20, 20 + maxlag):
        result["QSOgrant" + str(index + 1)] = round(
            weight * rightpart[index - 20] / temp, 6
        )

    return result


def GetQuasiSequenceOrder(
    ProteinSequence: str, maxlag: int = 30, weight: float = 0.1
) -> Dict[Any, Any]:
    """
    Compute quasi-sequence-order descriptors for a given protein.

    See [1]_ for details.

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence
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
    >>> from propy.GetProteinFromUniprot import GetProteinSequence
    >>> protein = GetProteinSequence(ProteinID="Q9NQ39")
    >>> result = GetQuasiSequenceOrder(protein)
    """
    result: Dict[Any, Any] = {}
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
    ProteinSequence: str,
    maxlag: int = 30,
    weight: float = 0.1,
    distancematrix: Dict[Any, Any] = None,
) -> Dict[Any, Any]:
    """
    Compute quasi-sequence-order descriptors for a given protein.

    See [1]_ for details.

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence
    maxlag : int, optional (default: 30)
        the maximum lag and the length of the protein should be larger than
        maxlag
    weight : float, optional (default: 0.1)
        a weight factor. Please see reference 1 for its choice.
    distancematrix : Dict[Any, Any]
        contains 400 distance values

    Returns
    -------
    result : Dict[Any, Any]
        contains all quasi-sequence-order descriptors

    Examples
    --------
    >>> from propy.GetProteinFromUniprot import GetProteinSequence
    >>> protein = GetProteinSequence(ProteinID="Q9NQ39")
    >>> result = GetQuasiSequenceOrderp(protein)
    """
    if distancematrix is None:
        distancematrix = {}
    result: Dict[Any, Any] = {}
    result.update(
        GetQuasiSequenceOrder1(ProteinSequence, maxlag, weight, distancematrix)
    )
    result.update(
        GetQuasiSequenceOrder2(ProteinSequence, maxlag, weight, distancematrix)
    )
    return result
