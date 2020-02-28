# -*- coding: utf-8 -*-
"""
This module is used for computing the Autocorrelation descriptors based
different properties of AADs. You can also input your properties of AADs, then
it can help you to compute Autocorrelation descriptors based on the property of
AADs. Currently, you can get 720 descriptors for a given protein sequence based
on our provided physicochemical properties of AADs.

References
----------
.. [1] http://www.genome.ad.jp/dbget/aaindex.html

.. [2] Feng, Z.P. and Zhang, C.T. (2000) Prediction of membrane protein types based on
   the hydrophobic index of amino acids. J Protein Chem, 19, 269-275.

.. [3] Horne, D.S. (1988) Prediction of protein helix content from an autocorrelation
   analysis of sequence hydrophobicities. Biopolymers, 27, 451-477.

.. [4] Sokal, R.R. and Thomson, B.A. (2006) Population structure inferred by local
   spatial autocorrelation: an Usage from an Amerindian tribal population. Am J
   Phys Anthropol, 129, 121-131.

Authors: Dongsheng Cao and Yizeng Liang.
Date: 2010.11.22
Email: oriental-cds@163.com
"""

# Core Library
import json
import math
from typing import Any, Dict, List

# Third party
from pkg_resources import resource_filename

AALetter: List[str] = list("ARNDCQEGHILKMFPSTWYV")

filepath = resource_filename(__name__, "data/hydrophobicity-autocorrelation.json")
with open(filepath, "r") as f:
    _Hydrophobicity: Dict[str, float] = json.load(f)

filepath = resource_filename(__name__, "data/AvFlexibility.json")
with open(filepath, "r") as f:
    _AvFlexibility: Dict[str, float] = json.load(f)

filepath = resource_filename(__name__, "data/Polarizability.json")
with open(filepath, "r") as f:
    _Polarizability: Dict[str, float] = json.load(f)

filepath = resource_filename(__name__, "data/FreeEnergy.json")
with open(filepath, "r") as f:
    _FreeEnergy: Dict[str, float] = json.load(f)

filepath = resource_filename(__name__, "data/ResidueASA.json")
with open(filepath, "r") as f:
    _ResidueASA: Dict[str, float] = json.load(f)


filepath = resource_filename(__name__, "data/ResidueVol.json")
with open(filepath, "r") as f:
    _ResidueVol: Dict[str, float] = json.load(f)


filepath = resource_filename(__name__, "data/Steric.json")
with open(filepath, "r") as f:
    _Steric: Dict[str, float] = json.load(f)

filepath = resource_filename(__name__, "data/Mutability.json")
with open(filepath, "r") as f:
    _Mutability: Dict[str, float] = json.load(f)


# Properties of AADs to compute the descriptors of protein sequence can
# continually be added.


_AAProperty = (
    _Hydrophobicity,
    _AvFlexibility,
    _Polarizability,
    _FreeEnergy,
    _ResidueASA,
    _ResidueVol,
    _Steric,
    _Mutability,
)

_AAPropertyName = (
    "_Hydrophobicity",
    "_AvFlexibility",
    "_Polarizability",
    "_FreeEnergy",
    "_ResidueASA",
    "_ResidueVol",
    "_Steric",
    "_Mutability",
)


def _mean(listvalue):
    """The mean value of the list data."""
    return sum(listvalue) / len(listvalue)


def _std(listvalue, ddof=1):
    """The standard deviation of the list data."""
    mean = _mean(listvalue)
    temp = [math.pow(i - mean, 2) for i in listvalue]
    res = math.sqrt(sum(temp) / (len(listvalue) - ddof))
    return res


def NormalizeEachAAP(AAP) -> Dict[Any, Any]:
    """
    Centralizes and standardizes all amino acid indices before the calculation

    Parameters
    ----------
    AAP : Dict
        contains the properties of 20 amino acids

    Returns
    -------
    result : Dict
        contains the normalized properties of 20 amino acids
    """
    if len(list(AAP.values())) != 20:
        print("You can not input the correct number of properities of Amino acids!")
    else:
        Result: Dict[Any, Any] = {}
        for i, j in list(AAP.items()):
            Result[i] = (j - _mean(list(AAP.values()))) / _std(
                list(AAP.values()), ddof=0
            )

    return Result


def CalculateEachNormalizedMoreauBrotoAuto(
    ProteinSequence: str, AAP: Dict[Any, Any], AAPName: str
) -> Dict[Any, Any]:
    """
    Compute MoreauBrotoAuto descriptors for different properties based on AADs.

    Examples
    --------
    >>> result = CalculateEachNormalizedMoreauBrotoAuto(protein,AAP,AAPName)

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence
    AAP : Dict[Any, Any]
        contains the properties of 20 amino acids (e.g., _AvFlexibility).
    AAPName : str
        used for indicating the property (e.g., '_AvFlexibility').

    Returns
    -------
    result contains 30 Normalized Moreau-Broto autocorrelation descriptors
    based on the given property.
    """
    AAPdic = NormalizeEachAAP(AAP)

    Result = {}
    for i in range(1, 31):
        temp = 0
        for j in range(len(ProteinSequence) - i):
            temp = temp + AAPdic[ProteinSequence[j]] * AAPdic[ProteinSequence[j + 1]]
        if len(ProteinSequence) - i == 0:
            Result["MoreauBrotoAuto" + AAPName + str(i)] = round(
                temp / (len(ProteinSequence)), 3
            )
        else:
            Result["MoreauBrotoAuto" + AAPName + str(i)] = round(
                temp / (len(ProteinSequence) - i), 3
            )

    return Result


def CalculateEachMoranAuto(
    ProteinSequence: str, AAP: Dict[Any, Any], AAPName: str
) -> Dict[Any, Any]:
    """
    Compute MoranAuto descriptors for different properties based on AADs.

    Examples
    --------
    >>> result = CalculateEachMoranAuto(protein,AAP,AAPName)

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence.
    AAP : Dict[Any, Any]
        contains the properties of 20 amino acids (e.g., _AvFlexibility).
    AAPName : str
        used for indicating the property (e.g., '_AvFlexibility').

    Returns
    -------
    result contains 30 Moran autocorrelation descriptors based on the given
    property.
    """
    AAPdic = NormalizeEachAAP(AAP)

    cds = 0
    for char in AALetter:
        cds = cds + (ProteinSequence.count(char)) * (AAPdic[char])
    Pmean = cds / len(ProteinSequence)

    cc = []
    for i in ProteinSequence:
        cc.append(AAPdic[i])

    K = (_std(cc, ddof=0)) ** 2

    Result = {}
    for i in range(1, 31):
        temp = 0
        for j in range(len(ProteinSequence) - i):
            temp = temp + (AAPdic[ProteinSequence[j]] - Pmean) * (
                AAPdic[ProteinSequence[j + i]] - Pmean
            )
        if len(ProteinSequence) - i == 0:
            Result["MoranAuto" + AAPName + str(i)] = round(
                temp / (len(ProteinSequence)) / K, 3
            )
        else:
            Result["MoranAuto" + AAPName + str(i)] = round(
                temp / (len(ProteinSequence) - i) / K, 3
            )

    return Result


def CalculateEachGearyAuto(
    ProteinSequence, AAP: Dict[Any, Any], AAPName
) -> Dict[Any, Any]:
    """
    Compute GearyAuto descriptors for different properties based on AADs.

    Examples
    --------
    >>> result = CalculateEachGearyAuto(protein,AAP,AAPName)

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence.
    AAP : Dict[Any, Any]
        contains the properties of 20 amino acids (e.g., _AvFlexibility).
    AAPName : str
        used for indicating the property (e.g., '_AvFlexibility').

    Returns
    -------
    result : Dict[Any, Any]
        contains 30 Geary autocorrelation descriptors based on the given
        property.
    """
    AAPdic = NormalizeEachAAP(AAP)

    cc = []
    for i in ProteinSequence:
        cc.append(AAPdic[i])

    K = ((_std(cc)) ** 2) * len(ProteinSequence) / (len(ProteinSequence) - 1)
    Result = {}
    for i in range(1, 31):
        temp = 0
        for j in range(len(ProteinSequence) - i):
            temp = (
                temp
                + (AAPdic[ProteinSequence[j]] - AAPdic[ProteinSequence[j + i]]) ** 2
            )
        if len(ProteinSequence) - i == 0:
            Result["GearyAuto" + AAPName + str(i)] = round(
                temp / (2 * (len(ProteinSequence))) / K, 3
            )
        else:
            Result["GearyAuto" + AAPName + str(i)] = round(
                temp / (2 * (len(ProteinSequence) - i)) / K, 3
            )
    return Result


def CalculateNormalizedMoreauBrotoAuto(
    ProteinSequence, AAProperty, AAPropertyName
) -> Dict[Any, Any]:
    """
    A method used for computing MoreauBrotoAuto for all properties.

    Examples
    --------
    >>> result = CalculateNormalizedMoreauBrotoAuto(protein,AAP,AAPName)

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence
    AAProperty : a list or tuple form
        contains the properties of 20 amino acids (e.g., _AAProperty).
    AAPName : a list or tuple form
        used for indicating the property (e.g., '_AAPropertyName').

    Returns
    -------
    result contains 30*p Normalized Moreau-Broto autocorrelation descriptors
    based on the given properties.
    """
    Result = {}
    for i in range(len(AAProperty)):
        Result[AAPropertyName[i]] = CalculateEachNormalizedMoreauBrotoAuto(
            ProteinSequence, AAProperty[i], AAPropertyName[i]
        )

    return Result


def CalculateMoranAuto(ProteinSequence, AAProperty, AAPropertyName) -> Dict[Any, Any]:
    """
    A method used for computing MoranAuto for all properties

    Examples
    --------
    >>> result = CalculateMoranAuto(protein,AAP,AAPName)

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence
    AAProperty : list or tuple form
        contains the properties of 20 amino acids (e.g., _AAProperty).
    AAPName : list or tuple form
        used for indicating the property (e.g., '_AAPropertyName').

    Returns
    -------
    result contains 30*p Moran autocorrelation descriptors based on the given
    properties.
    """
    Result = {}
    for i in range(len(AAProperty)):
        Result[AAPropertyName[i]] = CalculateEachMoranAuto(
            ProteinSequence, AAProperty[i], AAPropertyName[i]
        )

    return Result


def CalculateGearyAuto(ProteinSequence, AAProperty, AAPropertyName) -> Dict[Any, Any]:
    """
    A method used for computing GearyAuto for all properties.

    Examples
    --------
    >>> result = CalculateGearyAuto(protein, AAP, AAPName)

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    AAProperty : list or tuple form
        contains the properties of 20 amino acids (e.g., _AAProperty).
    AAPName : list or tuple form
        used for indicating the property (e.g., '_AAPropertyName').

    Returns
    -------
    result contains 30*p Geary autocorrelation descriptors based on the given
    properties.
    """
    Result = {}
    for i in range(len(AAProperty)):
        Result[AAPropertyName[i]] = CalculateEachGearyAuto(
            ProteinSequence, AAProperty[i], AAPropertyName[i]
        )

    return Result


# NormalizedMoreauBorto #######################################################
def CalculateNormalizedMoreauBrotoAutoHydrophobicity(ProteinSequence) -> Dict[Any, Any]:
    """
    Calculte the NormalizedMoreauBorto Autocorrelation descriptors based on
    hydrophobicity.

    Examples
    --------
    >>> result = CalculateNormalizedMoreauBrotoAutoHydrophobicity(protein)

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    Returns
    -------
    result contains 30 Normalized Moreau-Broto Autocorrelation descriptors
    based on Hydrophobicity.
    """
    result = CalculateEachNormalizedMoreauBrotoAuto(
        ProteinSequence, _Hydrophobicity, "_Hydrophobicity"
    )
    return result


def CalculateNormalizedMoreauBrotoAutoAvFlexibility(ProteinSequence) -> Dict[Any, Any]:
    """
    Calculte the NormalizedMoreauBorto Autocorrelation descriptors based on
    AvFlexibility.

    Examples
    --------
    >>> result = CalculateNormalizedMoreauBrotoAutoAvFlexibility(protein)

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    Returns
    -------
    result contains 30 Normalized Moreau-Broto Autocorrelation descriptors
    based on AvFlexibility.
    """
    result = CalculateEachNormalizedMoreauBrotoAuto(
        ProteinSequence, _AvFlexibility, "_AvFlexibility"
    )
    return result


def CalculateNormalizedMoreauBrotoAutoPolarizability(
    ProteinSequence: str,
) -> Dict[Any, Any]:
    """
    Calculte the NormalizedMoreauBorto Autocorrelation descriptors based on
    Polarizability.

    Examples
    --------
    >>> result = CalculateNormalizedMoreauBrotoAutoPolarizability(protein)

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    Returns
    -------
    result contains 30 Normalized Moreau-Broto Autocorrelation descriptors
    based on Polarizability.
    """
    result = CalculateEachNormalizedMoreauBrotoAuto(
        ProteinSequence, _Polarizability, "_Polarizability"
    )
    return result


def CalculateNormalizedMoreauBrotoAutoFreeEnergy(
    ProteinSequence: str,
) -> Dict[Any, Any]:
    """
    Calculte the NormalizedMoreauBorto Autocorrelation descriptors based on
    FreeEnergy.

    Examples
    --------
    >>> result = CalculateNormalizedMoreauBrotoAutoFreeEnergy(protein)

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    Returns
    -------
    result contains 30 Normalized Moreau-Broto Autocorrelation descriptors
    based on FreeEnergy.
    """
    result = CalculateEachNormalizedMoreauBrotoAuto(
        ProteinSequence, _FreeEnergy, "_FreeEnergy"
    )
    return result


def CalculateNormalizedMoreauBrotoAutoResidueASA(
    ProteinSequence: str,
) -> Dict[Any, Any]:
    """
    Calculte the NormalizedMoreauBorto Autocorrelation descriptors based on
    ResidueASA.

    Examples
    --------
    >>> result = CalculateNormalizedMoreauBrotoAutoResidueASA(protein)

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    Returns
    -------
    result contains 30 Normalized Moreau-Broto Autocorrelation descriptors
    based on ResidueASA.
    """
    result = CalculateEachNormalizedMoreauBrotoAuto(
        ProteinSequence, _ResidueASA, "_ResidueASA"
    )
    return result


def CalculateNormalizedMoreauBrotoAutoResidueVol(ProteinSequence) -> Dict[Any, Any]:
    """
    Calculte the NormalizedMoreauBorto Autocorrelation descriptors based on
    ResidueVol.

    Examples
    --------
    >>> result = CalculateNormalizedMoreauBrotoAutoResidueVol(protein)

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    Returns
    -------
    result contains 30 Normalized Moreau-Broto Autocorrelation descriptors
    based on ResidueVol.
    """
    result = CalculateEachNormalizedMoreauBrotoAuto(
        ProteinSequence, _ResidueVol, "_ResidueVol"
    )
    return result


def CalculateNormalizedMoreauBrotoAutoSteric(ProteinSequence) -> Dict[Any, Any]:
    """
    Calculte the NormalizedMoreauBorto Autocorrelation descriptors based on Steric.

    Examples
    --------
    >>> result = CalculateNormalizedMoreauBrotoAutoSteric(protein)

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    Returns
    -------
    result contains 30 Normalized Moreau-Broto Autocorrelation descriptors
    based on Steric.
    """
    result = CalculateEachNormalizedMoreauBrotoAuto(ProteinSequence, _Steric, "_Steric")
    return result


def CalculateNormalizedMoreauBrotoAutoMutability(ProteinSequence) -> Dict[Any, Any]:
    """
    Calculte the NormalizedMoreauBorto Autocorrelation descriptors based on Mutability.

    Examples
    --------
    >>> result = CalculateNormalizedMoreauBrotoAutoMutability(protein)

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    Returns
    -------
    result contains 30 Normalized Moreau-Broto Autocorrelation descriptors
    based on Mutability.
    """
    result = CalculateEachNormalizedMoreauBrotoAuto(
        ProteinSequence, _Mutability, "_Mutability"
    )
    return result


# MoranAuto ###################################################################
def CalculateMoranAutoHydrophobicity(ProteinSequence) -> Dict[Any, Any]:
    """
    Calculte the MoranAuto Autocorrelation descriptors based on hydrophobicity.

    Examples
    --------
    >>> result = CalculateMoranAutoHydrophobicity(protein)

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    Returns
    -------
    result contains 30 Moran Autocorrelation descriptors based on
    hydrophobicity.
    """
    result = CalculateEachMoranAuto(ProteinSequence, _Hydrophobicity, "_Hydrophobicity")
    return result


def CalculateMoranAutoAvFlexibility(ProteinSequence: str) -> Dict[Any, Any]:
    """
    Calculte the MoranAuto Autocorrelation descriptors based on AvFlexibility.

    Examples
    --------
    >>> result = CalculateMoranAutoAvFlexibility(protein)

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    Returns
    -------
    result contains 30 Moran Autocorrelation descriptors based on
    AvFlexibility.
    """
    result = CalculateEachMoranAuto(ProteinSequence, _AvFlexibility, "_AvFlexibility")
    return result


def CalculateMoranAutoPolarizability(ProteinSequence: str) -> Dict[Any, Any]:
    """
    Calculte the MoranAuto Autocorrelation descriptors based on Polarizability.

    Examples
    --------
    >>> result = CalculateMoranAutoPolarizability(protein)

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    Returns
    -------
    result contains 30 Moran Autocorrelation descriptors based on
    Polarizability.
    """
    result = CalculateEachMoranAuto(ProteinSequence, _Polarizability, "_Polarizability")
    return result


def CalculateMoranAutoFreeEnergy(ProteinSequence: str) -> Dict[Any, Any]:
    """
    Calculte the MoranAuto Autocorrelation descriptors based on

    FreeEnergy.

    Examples
    --------
    >>> result = CalculateMoranAutoFreeEnergy(protein)

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    Returns
    -------
    result contains 30 Moran Autocorrelation descriptors based on FreeEnergy.
    """
    result = CalculateEachMoranAuto(ProteinSequence, _FreeEnergy, "_FreeEnergy")
    return result


def CalculateMoranAutoResidueASA(ProteinSequence: str) -> Dict[Any, Any]:
    """
    Calculte the MoranAuto Autocorrelation descriptors based on ResidueASA.

    Examples
    --------
    >>> result = CalculateMoranAutoResidueASA(protein)

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    Returns
    -------
    result contains 30 Moran Autocorrelation descriptors based on ResidueASA.
    """
    result = CalculateEachMoranAuto(ProteinSequence, _ResidueASA, "_ResidueASA")
    return result


def CalculateMoranAutoResidueVol(ProteinSequence: str) -> Dict[Any, Any]:
    """
    Calculte the MoranAuto Autocorrelation descriptors based on ResidueVol.

    Examples
    --------
    >>> result = CalculateMoranAutoResidueVol(protein)

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    Returns
    -------
    result contains 30 Moran Autocorrelation descriptors based on ResidueVol.
    """
    result = CalculateEachMoranAuto(ProteinSequence, _ResidueVol, "_ResidueVol")
    return result


def CalculateMoranAutoSteric(ProteinSequence: str) -> Dict[Any, Any]:
    """
    Calculte the MoranAuto Autocorrelation descriptors based on AutoSteric.

    Examples
    --------
    >>> result = CalculateMoranAutoSteric(protein)

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    Returns
    -------
    result contains 30 Moran Autocorrelation descriptors based on AutoSteric.
    """
    result = CalculateEachMoranAuto(ProteinSequence, _Steric, "_Steric")
    return result


def CalculateMoranAutoMutability(ProteinSequence: str) -> Dict[Any, Any]:
    """
    Calculte the MoranAuto Autocorrelation descriptors based on Mutability.

    Examples
    --------
    >>> result = CalculateMoranAutoMutability(protein)

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    Returns
    -------
    result contains 30 Moran Autocorrelation descriptors based on Mutability.
    """
    result = CalculateEachMoranAuto(ProteinSequence, _Mutability, "_Mutability")
    return result


# GearyAuto####################################################################
def CalculateGearyAutoHydrophobicity(ProteinSequence: str) -> Dict[Any, Any]:
    """
    Calculte the GearyAuto Autocorrelation descriptors based on hydrophobicity.

    Examples
    --------
    >>> result = CalculateGearyAutoHydrophobicity(protein)

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    Returns
    -------
    result contains 30 Geary Autocorrelation descriptors based on
    hydrophobicity.
    """
    result = CalculateEachGearyAuto(ProteinSequence, _Hydrophobicity, "_Hydrophobicity")
    return result


def CalculateGearyAutoAvFlexibility(ProteinSequence):
    """
    Calculte the GearyAuto Autocorrelation descriptors based on AvFlexibility.

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence.

    Returns
    -------
    result : Dict
        contains 30 Geary Autocorrelation descriptors based on AvFlexibility.

    Examples
    --------
    >>> result = CalculateGearyAutoAvFlexibility(protein)
    """
    result = CalculateEachGearyAuto(ProteinSequence, _AvFlexibility, "_AvFlexibility")
    return result


def CalculateGearyAutoPolarizability(ProteinSequence):
    """
    Calculte the GearyAuto Autocorrelation descriptors based on Polarizability.

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence.

    Returns
    -------
    result : Dict
        contains 30 Geary Autocorrelation descriptors based on Polarizability.

    Examples
    --------
    >>> result = CalculateGearyAutoPolarizability(protein)
    """
    result = CalculateEachGearyAuto(ProteinSequence, _Polarizability, "_Polarizability")
    return result


def CalculateGearyAutoFreeEnergy(ProteinSequence: str) -> Dict[Any, Any]:
    """
    Calculte the GearyAuto Autocorrelation descriptors based on FreeEnergy.

    Examples
    --------
    >>> result = CalculateGearyAutoFreeEnergy(protein)

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    Returns
    -------
    result contains 30 Geary Autocorrelation descriptors based on FreeEnergy.
    """
    result = CalculateEachGearyAuto(ProteinSequence, _FreeEnergy, "_FreeEnergy")
    return result


def CalculateGearyAutoResidueASA(ProteinSequence: str) -> Dict[Any, Any]:
    """
    Calculte the GearyAuto Autocorrelation descriptors based on ResidueASA.

    Examples
    --------
    >>> result = CalculateGearyAutoResidueASA(protein)

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    Returns
    -------
    result contains 30 Geary Autocorrelation descriptors based on ResidueASA.
    """
    result = CalculateEachGearyAuto(ProteinSequence, _ResidueASA, "_ResidueASA")
    return result


def CalculateGearyAutoResidueVol(ProteinSequence: str) -> Dict[Any, Any]:
    """
    Calculte the GearyAuto Autocorrelation descriptors based on ResidueVol.

    Examples
    --------
    >>> result = CalculateGearyAutoResidueVol(protein)

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    Returns
    -------
    result contains 30 Geary Autocorrelation descriptors based on ResidueVol.
    """
    result = CalculateEachGearyAuto(ProteinSequence, _ResidueVol, "_ResidueVol")
    return result


def CalculateGearyAutoSteric(ProteinSequence: str) -> Dict[Any, Any]:
    """
    Calculte the GearyAuto Autocorrelation descriptors based on Steric.

    Examples
    --------
    >>> result = CalculateGearyAutoSteric(protein)

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    Returns
    -------
    result contains 30 Geary Autocorrelation descriptors based on Steric.
    """
    result = CalculateEachGearyAuto(ProteinSequence, _Steric, "_Steric")
    return result


def CalculateGearyAutoMutability(ProteinSequence: str) -> Dict[Any, Any]:
    """
    Calculte the GearyAuto Autocorrelation descriptors based on Mutability.

    Examples
    --------
    >>> result = CalculateGearyAutoMutability(protein)

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    Returns
    -------
    result contains 30 Geary Autocorrelation descriptors based on Mutability.
    """
    result = CalculateEachGearyAuto(ProteinSequence, _Mutability, "_Mutability")
    return result


def CalculateNormalizedMoreauBrotoAutoTotal(ProteinSequence: str) -> Dict[Any, Any]:
    """
    A method used for computing normalized Moreau Broto autocorrelation
    descriptors based on 8 proterties of AADs.

    Examples
    --------
    >>> result = CalculateNormalizedMoreauBrotoAutoTotal(protein)

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    Returns
    -------
    result contains 30*8=240 normalized Moreau Broto autocorrelation
    descriptors based on the given properties(i.e., _AAPropert).
    """
    result: Dict[Any, Any] = {}
    result.update(CalculateNormalizedMoreauBrotoAutoHydrophobicity(ProteinSequence))
    result.update(CalculateNormalizedMoreauBrotoAutoAvFlexibility(ProteinSequence))
    result.update(CalculateNormalizedMoreauBrotoAutoPolarizability(ProteinSequence))
    result.update(CalculateNormalizedMoreauBrotoAutoFreeEnergy(ProteinSequence))
    result.update(CalculateNormalizedMoreauBrotoAutoResidueASA(ProteinSequence))
    result.update(CalculateNormalizedMoreauBrotoAutoResidueVol(ProteinSequence))
    result.update(CalculateNormalizedMoreauBrotoAutoSteric(ProteinSequence))
    result.update(CalculateNormalizedMoreauBrotoAutoMutability(ProteinSequence))
    return result


def CalculateMoranAutoTotal(ProteinSequence: str) -> Dict[Any, Any]:
    """
    A method used for computing Moran autocorrelation descriptors based on 8
    properties of AADs.

    Examples
    --------
    >>> result = CalculateMoranAutoTotal(protein)

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    Returns
    -------
    result contains 30*8=240 Moran autocorrelation descriptors based on the
    given properties(i.e., _AAPropert).
    """
    result: Dict[Any, Any] = {}
    result.update(CalculateMoranAutoHydrophobicity(ProteinSequence))
    result.update(CalculateMoranAutoAvFlexibility(ProteinSequence))
    result.update(CalculateMoranAutoPolarizability(ProteinSequence))
    result.update(CalculateMoranAutoFreeEnergy(ProteinSequence))
    result.update(CalculateMoranAutoResidueASA(ProteinSequence))
    result.update(CalculateMoranAutoResidueVol(ProteinSequence))
    result.update(CalculateMoranAutoSteric(ProteinSequence))
    result.update(CalculateMoranAutoMutability(ProteinSequence))
    return result


def CalculateGearyAutoTotal(ProteinSequence: str) -> Dict[Any, Any]:
    """
    A method used for computing Geary autocorrelation descriptors based on 8
    properties of AADs.

    Examples
    --------
    >>> result = CalculateGearyAutoTotal(protein)

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    Returns
    -------
    result contains 30*8=240 Geary autocorrelation descriptors based on the
    given properties(i.e., _AAPropert).
    """
    result: Dict[Any, Any] = {}
    result.update(CalculateGearyAutoHydrophobicity(ProteinSequence))
    result.update(CalculateGearyAutoAvFlexibility(ProteinSequence))
    result.update(CalculateGearyAutoPolarizability(ProteinSequence))
    result.update(CalculateGearyAutoFreeEnergy(ProteinSequence))
    result.update(CalculateGearyAutoResidueASA(ProteinSequence))
    result.update(CalculateGearyAutoResidueVol(ProteinSequence))
    result.update(CalculateGearyAutoSteric(ProteinSequence))
    result.update(CalculateGearyAutoMutability(ProteinSequence))
    return result


def CalculateAutoTotal(ProteinSequence: str) -> Dict[Any, Any]:
    """
    A method used for computing all autocorrelation descriptors based on 8
    properties of AADs.

    Examples
    --------
    >>> result = CalculateGearyAutoTotal(protein)

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    Returns
    -------
    result contains 30*8*3=720 normalized Moreau Broto, Moran, and Geary
    autocorrelation descriptors based on the given properties(i.e.,
    _AAPropert).
    """
    result: Dict[Any, Any] = {}
    result.update(CalculateNormalizedMoreauBrotoAutoTotal(ProteinSequence))
    result.update(CalculateMoranAutoTotal(ProteinSequence))
    result.update(CalculateGearyAutoTotal(ProteinSequence))
    return result
