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
Instead of using the conventional 20-D amino acid composition to represent the
sample of a protein, Prof. Kuo-Chen Chou proposed the pseudo amino acid (PseAA)
composition in order for inluding the sequence-order information. Based on the
concept of Chou's pseudo amino acid composition, the server PseAA was designed
in a flexible way, allowing users to generate various kinds of pseudo amino
acid composition for a given protein sequence by selecting different parameters
and their combinations. This module aims at computing two types of PseAA
descriptors: Type I and Type II.

References
----------
.. [1] Kuo-Chen Chou. Prediction of Protein Cellular Attributes Using
       Pseudo-Amino Acid Composition. PROTEINS: Structure, Function, and
       Genetics, 2001, 43: 246-255.

.. [2] http://www.csbio.sjtu.edu.cn/bioinf/PseAAC/
.. [3] http://www.csbio.sjtu.edu.cn/bioinf/PseAAC/type2.htm
.. [4] Kuo-Chen Chou. Using amphiphilic pseudo amino acid composition to
       predict enzyme subfamily classes. Bioinformatics, 2005, 21, 10-19.


The hydrophobicity values are from JACS, 1962, 84: 4240-4246. (C. Tanford).

The hydrophilicity values are from PNAS, 1981, 78:3824-3828 (T.P.Hopp & K.R.Woods).

The side-chain mass for each of the 20 amino acids.

CRC Handbook of Chemistry and Physics, 66th ed., CRC Press, Boca Raton, Florida (1985).

R.M.C. Dawson, D.C. Elliott, W.H. Elliott, K.M. Jones,
Data for Biochemical Research 3rd ed.,

Clarendon Press Oxford (1986).
"""

# Core Library
import json
import math
from typing import Any, Dict

# Third party
from pkg_resources import resource_filename

# First party
from propy import AALetter

with open(resource_filename(__name__, "data/hydrophobicity.json"), "r") as f:
    _Hydrophobicity: Dict[str, float] = json.load(f)

with open(resource_filename(__name__, "data/hydrophilicity.json"), "r") as f:
    _hydrophilicity: Dict[str, float] = json.load(f)

with open(resource_filename(__name__, "data/residuemass.json"), "r") as f:
    _residuemass: Dict[str, float] = json.load(f)


with open(resource_filename(__name__, "data/pK1.json"), "r") as f:
    _pK1: Dict[str, float] = json.load(f)

with open(resource_filename(__name__, "data/pK2.json"), "r") as f:
    _pK2: Dict[str, float] = json.load(f)

with open(resource_filename(__name__, "data/pI.json"), "r") as f:
    _pI: Dict[str, float] = json.load(f)


def _mean(listvalue):
    """
    The mean value of the list data.

    Examples
    --------
    >>> _mean(listvalue=[1, 2, 3])
    2.0
    """
    return sum(listvalue) / len(listvalue)


def _std(listvalue, ddof=1):
    """
    The standard deviation of the list data.

    Examples
    --------
    >>> _std(listvalue=[1, 2, 3])
    1.0
    """
    mean = _mean(listvalue)
    temp = [math.pow(i - mean, 2) for i in listvalue]
    res = math.sqrt(sum(temp) / (len(listvalue) - ddof))
    return res


def NormalizeEachAAP(AAP):
    """
    All of the amino acid indices are centralized and standardized before the
    calculation.

    Parameters
    ----------
    AAP is a dict form containing the properties of 20 amino acids.

    Returns
    -------
    result is the a dict form containing the normalized properties of 20 amino
    acids.

    Examples
    --------
    >>> result = NormalizeEachAAP(AAP=_Hydrophobicity)
    """
    if len(list(AAP.values())) != 20:
        print("You can not input the correct number of properities of Amino acids!")
    else:
        result = {}
        for i, j in list(AAP.items()):
            result[i] = (j - _mean(list(AAP.values()))) / _std(
                list(AAP.values()), ddof=0
            )

    return result


# Type I descriptors###########################################################
# Pseudo-Amino Acid Composition descriptors####################################
def _GetCorrelationFunction(
    Ri="S", Rj="D", AAP=(_Hydrophobicity, _hydrophilicity, _residuemass)
):
    """
    Computing the correlation between two given amino acids using the above
    three properties.

    Parameters
    ----------
    Ri and Rj are the amino acids, respectively.

    Returns
    -------
    result is the correlation value between two amino acids.

    Examples
    --------
    >>> result = _GetCorrelationFunction(Ri="S", Rj="D")
    """
    Hydrophobicity = NormalizeEachAAP(AAP[0])
    hydrophilicity = NormalizeEachAAP(AAP[1])
    residuemass = NormalizeEachAAP(AAP[2])
    theta1 = math.pow(Hydrophobicity[Ri] - Hydrophobicity[Rj], 2)
    theta2 = math.pow(hydrophilicity[Ri] - hydrophilicity[Rj], 2)
    theta3 = math.pow(residuemass[Ri] - residuemass[Rj], 2)
    theta = round((theta1 + theta2 + theta3) / 3.0, 3)
    return theta


def _GetSequenceOrderCorrelationFactor(ProteinSequence: str, k: int = 1) -> float:
    """
    Computing the Sequence order correlation factor with gap equal to k based
    on [_Hydrophobicity, _hydrophilicity, _residuemass].

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence.
    k : int
        is the gap.

    Returns
    -------
    result : float
        the correlation factor value with the gap equal to k

    Examples
    --------
    >>> from propy.GetProteinFromUniprot import GetProteinSequence
    >>> protein = GetProteinSequence(ProteinID="Q9NQ39")
    >>> result = _GetSequenceOrderCorrelationFactor(protein)
    """
    LengthSequence = len(ProteinSequence)
    res = []
    for i in range(LengthSequence - k):
        AA1 = ProteinSequence[i]
        AA2 = ProteinSequence[i + k]
        res.append(_GetCorrelationFunction(AA1, AA2))
    result = round(sum(res) / (LengthSequence - k), 3)
    return result


def GetAAComposition(ProteinSequence: str) -> Dict[Any, Any]:
    """
    Calculate the composition of Amino acids for a given protein sequence.

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    Returns
    -------
    result is a dict form containing the composition of 20 amino acids.

    Examples
    --------
    >>> from propy.GetProteinFromUniprot import GetProteinSequence
    >>> protein = GetProteinSequence(ProteinID="Q9NQ39")
    >>> result = GetAAComposition(protein)
    """
    LengthSequence = len(ProteinSequence)
    Result = {}
    for i in AALetter:
        Result[i] = round(float(ProteinSequence.count(i)) / LengthSequence * 100, 3)
    return Result


def _GetPseudoAAC1(ProteinSequence, lamda=10, weight=0.05):
    """
    Computing the first 20 of type I pseudo-amino acid compostion descriptors based on

    [_Hydrophobicity, _hydrophilicity, _residuemass].
    """
    rightpart = 0.0
    for i in range(lamda):
        rightpart = rightpart + _GetSequenceOrderCorrelationFactor(
            ProteinSequence, k=i + 1
        )
    AAC = GetAAComposition(ProteinSequence)

    result = {}
    temp = 1 + weight * rightpart
    for index, char in enumerate(AALetter):
        result["PAAC" + str(index + 1)] = round(AAC[char] / temp, 3)

    return result


def _GetPseudoAAC2(ProteinSequence, lamda=10, weight=0.05):
    """
    Computing the last lamda of type I pseudo-amino acid compostion descriptors based on

    [_Hydrophobicity, _hydrophilicity, _residuemass].
    """
    rightpart = []
    for i in range(lamda):
        rightpart.append(_GetSequenceOrderCorrelationFactor(ProteinSequence, k=i + 1))

    result = {}
    temp = 1 + weight * sum(rightpart)
    for index in range(20, 20 + lamda):
        result["PAAC" + str(index + 1)] = round(
            weight * rightpart[index - 20] / temp * 100, 3
        )

    return result


def _GetPseudoAAC(
    ProteinSequence: str, lamda: int = 10, weight: float = 0.05
) -> Dict[Any, Any]:
    """
    Computing all of type I pseudo-amino acid compostion descriptors based on
    three given properties. Note that the number of PAAC strongly depends on
    the lamda value. if lamda = 20, we can obtain 20+20=40 PAAC descriptors.
    The size of these values depends on the choice of lamda and weight
    simultaneously.

    AAP=[_Hydrophobicity, _hydrophilicity, _residuemass]

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence
    lamda : int
        reflects the rank of correlation and is a non-Negative integer, such as
        15. Note that (1) lamda should NOT be larger than the length of input
        protein sequence; (2) lamda must be non-Negative integer, such as 0, 1,
        2, ...; (3) when lamda =0, the output of PseAA server is the 20-D amino
        acid composition.
    weight factor : float
        is designed for the users to put weight on the additional PseAA
        components with respect to the conventional AA components. The user can
        select any value within the region from 0.05 to 0.7 for the weight
        factor.

    Returns
    -------
    result : Dict[Any, Any]
        contains calculated 20+lamda PAAC descriptors

    Examples
    --------
    >>> from propy.GetProteinFromUniprot import GetProteinSequence
    >>> protein = GetProteinSequence(ProteinID="Q9NQ39")
    >>> result = _GetPseudoAAC(protein)
    """
    res: Dict[Any, Any] = {}
    res.update(_GetPseudoAAC1(ProteinSequence, lamda=lamda, weight=weight))
    res.update(_GetPseudoAAC2(ProteinSequence, lamda=lamda, weight=weight))
    return res


# Type II descriptors##########################################################
# Amphiphilic Pseudo-Amino Acid Composition descriptors########################
def _GetCorrelationFunctionForAPAAC(
    Ri="S", Rj="D", AAP=(_Hydrophobicity, _hydrophilicity)
):
    """
    Computing the correlation between two given amino acids using the above two
    properties for APAAC (type II PseAAC).

    Parameters
    ----------
    Ri and Rj are the amino acids, respectively.

    Returns
    -------
    result :
        the correlation value between two amino acids

    Examples
    --------
    >>> result = _GetCorrelationFunctionForAPAAC(Ri="S", Rj="D")
    """
    Hydrophobicity = NormalizeEachAAP(AAP[0])
    hydrophilicity = NormalizeEachAAP(AAP[1])
    theta1 = round(Hydrophobicity[Ri] * Hydrophobicity[Rj], 3)
    theta2 = round(hydrophilicity[Ri] * hydrophilicity[Rj], 3)

    return theta1, theta2


def GetSequenceOrderCorrelationFactorForAPAAC(ProteinSequence, k=1):
    """
    Computing the Sequence order correlation factor with gap equal to k based on

    [_Hydrophobicity, _hydrophilicity] for APAAC (type II PseAAC) .

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    k is the gap.

    Returns
    -------
    result is the correlation factor value with the gap equal to k.

    Examples
    --------
    >>> from propy.GetProteinFromUniprot import GetProteinSequence
    >>> protein = GetProteinSequence(ProteinID="Q9NQ39")
    >>> result = GetSequenceOrderCorrelationFactorForAPAAC(protein)
    """
    LengthSequence = len(ProteinSequence)
    resHydrophobicity = []
    reshydrophilicity = []
    for i in range(LengthSequence - k):
        AA1 = ProteinSequence[i]
        AA2 = ProteinSequence[i + k]
        temp = _GetCorrelationFunctionForAPAAC(AA1, AA2)
        resHydrophobicity.append(temp[0])
        reshydrophilicity.append(temp[1])
    result = []
    result.append(round(sum(resHydrophobicity) / (LengthSequence - k), 3))
    result.append(round(sum(reshydrophilicity) / (LengthSequence - k), 3))
    return result


def GetAPseudoAAC1(ProteinSequence, lamda=30, weight=0.5):
    """
    Computing the first 20 of type II pseudo-amino acid compostion descriptors based on

    [_Hydrophobicity, _hydrophilicity].
    """
    rightpart = 0.0
    for i in range(lamda):
        rightpart = rightpart + sum(
            GetSequenceOrderCorrelationFactorForAPAAC(ProteinSequence, k=i + 1)
        )
    AAC = GetAAComposition(ProteinSequence)

    result = {}
    temp = 1 + weight * rightpart
    for index, char in enumerate(AALetter):
        result["APAAC" + str(index + 1)] = round(AAC[char] / temp, 3)

    return result


def GetAPseudoAAC2(ProteinSequence, lamda=30, weight=0.5):
    """
    Computing the last lamda of type II pseudo-amino acid compostion descriptors

    based on (_Hydrophobicity, _hydrophilicity).
    """
    rightpart = []
    for i in range(lamda):
        temp = GetSequenceOrderCorrelationFactorForAPAAC(ProteinSequence, k=i + 1)
        rightpart.append(temp[0])
        rightpart.append(temp[1])

    result = {}
    temp = 1 + weight * sum(rightpart)
    for index in range(20, 20 + 2 * lamda):
        result["APAAC" + str(index + 1)] = round(
            weight * rightpart[index - 20] / temp * 100, 3
        )

    return result


def GetAPseudoAAC(ProteinSequence, lamda: int = 30, weight: float = 0.5):
    """
    Computing all of type II pseudo-amino acid compostion descriptors based on
    the given properties. Note that the number of PAAC strongly depends on the
    lamda value. if lamda = 20, we can obtain 20+20=40 PAAC descriptors. The
    size of these values depends on the choice of lamda and weight
    simultaneously.

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    lamda : int
        reflects the rank of correlation and is a non-Negative integer, such as
        15. Note that (1)lamda should NOT be larger than the length of input
        protein sequence; (2) lamda must be non-Negative integer, such as 0, 1,
        2, ...; (3) when lamda =0, the output of PseAA server is the 20-D amino
        acid composition.
    weight : float
        is designed for the users to put weight on the additional PseAA
        components with respect to the conventional AA components. The user can
        select any value within the region from 0.05 to 0.7 for the weight
        factor.

    Returns
    -------
    result : Dict[Any, Any]
        contains calculated 20+lamda PAAC descriptors

    Examples
    --------
    >>> from propy.GetProteinFromUniprot import GetProteinSequence
    >>> protein = GetProteinSequence(ProteinID="Q9NQ39")
    >>> result = GetAPseudoAAC(protein)
    """
    res: Dict[Any, Any] = {}
    res.update(GetAPseudoAAC1(ProteinSequence, lamda=lamda, weight=weight))
    res.update(GetAPseudoAAC2(ProteinSequence, lamda=lamda, weight=weight))
    return res


# Type I descriptors###########################################################
# Pseudo-Amino Acid Composition descriptors####################################
# based on different properties################################################
def GetCorrelationFunction(Ri="S", Rj="D", AAP=None):
    """
    Computing the correlation between two given amino acids using the given
    properties.

    Parameters
    ----------
    Ri : str
        amino acids
    Rj : str
        amino acids
    AAP : List[Any]
        contains the properties, each of which is a dict form.

    Returns
    -------
    result is the correlation value between two amino acids.

    Examples
    --------
    >>> GetCorrelationFunction(Ri="S", Rj="D", AAP=_Hydrophobicity)
    """
    if AAP is None:
        AAP = []
    NumAAP = len(AAP)
    theta = 0.0
    for i in range(NumAAP):
        temp = NormalizeEachAAP(AAP[i])
        theta = theta + math.pow(temp[Ri] - temp[Rj], 2)
    result = round(theta / NumAAP, 3)
    return result


def GetSequenceOrderCorrelationFactor(ProteinSequence, k: int = 1, AAP=None):
    """
    Computing the Sequence order correlation factor with gap equal to k based
    on the given properities.

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence
    k : int
        the gap.
    AAP : List[Any]
        contains the properties, each of which is a dict form.

    Returns
    -------
    result is the correlation factor value with the gap equal to k.

    Examples
    --------
    >>> from propy.GetProteinFromUniprot import GetProteinSequence
    >>> protein = GetProteinSequence(ProteinID="Q9NQ39")
    >>> result = GetSequenceOrderCorrelationFactor(protein)
    """
    if AAP is None:
        AAP = []
    LengthSequence = len(ProteinSequence)
    res = []
    for i in range(LengthSequence - k):
        AA1 = ProteinSequence[i]
        AA2 = ProteinSequence[i + k]
        res.append(GetCorrelationFunction(AA1, AA2, AAP))
    result = round(sum(res) / (LengthSequence - k), 3)
    return result


def GetPseudoAAC1(ProteinSequence, lamda=30, weight=0.05, AAP=None):
    """
    Computing the first 20 of type I pseudo-amino acid compostion descriptors
    based on the given properties.
    """
    if AAP is None:
        AAP = []
    rightpart = 0.0
    for i in range(lamda):
        rightpart = rightpart + GetSequenceOrderCorrelationFactor(
            ProteinSequence, i + 1, AAP
        )
    AAC = GetAAComposition(ProteinSequence)

    result = {}
    temp = 1 + weight * rightpart
    for index, char in enumerate(AALetter):
        result["PAAC" + str(index + 1)] = round(AAC[char] / temp, 3)

    return result


def GetPseudoAAC2(ProteinSequence, lamda: int = 30, weight: float = 0.05, AAP=None):
    """
    Compute the last lamda of type I pseudo-amino acid compostion descriptors
    based on the given properties.
    """
    if AAP is None:
        AAP = []
    rightpart = []
    for i in range(lamda):
        rightpart.append(GetSequenceOrderCorrelationFactor(ProteinSequence, i + 1, AAP))

    result = {}
    temp = 1 + weight * sum(rightpart)
    for index in range(20, 20 + lamda):
        result["PAAC" + str(index + 1)] = round(
            weight * rightpart[index - 20] / temp * 100, 3
        )

    return result


def GetPseudoAAC(ProteinSequence: str, lamda: int = 30, weight: float = 0.05, AAP=None):
    """
    Computing all of type I pseudo-amino acid compostion descriptors based on
    the given properties. Note that the number of PAAC strongly depends on the
    lamda value. if lamda = 20, we can obtain 20+20=40 PAAC descriptors. The
    size of these values depends on the choice of lamda and weight
    simultaneously. You must specify some properties into AAP.

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence
    lamda : int
        reflects the rank of correlation and is a non-Negative integer, such as
        15. Note that (1)lamda should NOT be larger than the length of input
        protein sequence; (2) lamda must be non-Negative integer, such as 0, 1,
        2, ...; (3) when lamda =0, the output of PseAA server is the 20-D amino
        acid composition.
    weight : float
        is designed for the users to put weight on the additional PseAA components
        with respect to the conventional AA components. The user can select any
        value within the region from 0.05 to 0.7 for the weight factor.
    AAP : List[Any]
        contains the properties, each of which is a dict form.

    Returns
    -------
    result is a dict form containing calculated 20+lamda PAAC descriptors.

    Examples
    --------
    >>> from propy.GetProteinFromUniprot import GetProteinSequence
    >>> protein = GetProteinSequence(ProteinID="Q9NQ39")
    >>> result = GetPseudoAAC(protein)
    """
    if AAP is None:
        AAP = []
    res: Dict[Any, Any] = {}
    res.update(GetPseudoAAC1(ProteinSequence, lamda, weight, AAP))
    res.update(GetPseudoAAC2(ProteinSequence, lamda, weight, AAP))
    return res
