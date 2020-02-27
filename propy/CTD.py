# -*- coding: utf-8 -*-
"""
This module is used for computing the composition, transition and distribution
descriptors based on the different properties of AADs. The AADs with the same
properties is marked as the same number. You can get 147 descriptors for a
given protein sequence.

References
----------
.. [1] Inna Dubchak, Ilya Muchink, Stephen R.Holbrook and Sung-Hou Kim.
       Prediction of protein folding class using global description of amino
       acid sequence. Proc.Natl. Acad.Sci.USA, 1995, 92, 8700-8704.

.. [2] Inna Dubchak, Ilya Muchink, Christopher Mayor, Igor Dralyuk and Sung-Hou
       Kim. Recognition of a Protein Fold in the Context of the SCOP
       classification. Proteins: Structure, Function and
       Genetics,1999,35,401-407.

Authors: Dongsheng Cao and Yizeng Liang.
Date: 2010.11.22
Email: oriental-cds@163.com
"""

# Core Library
import copy
import math
from typing import Any, Dict

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

_Hydrophobicity = {"1": "RKEDQN", "2": "GASTPHY", "3": "CLVIMFW"}
# '1'stand for Polar; '2'stand for Neutral, '3' stand for Hydrophobicity

_NormalizedVDWV = {"1": "GASTPD", "2": "NVEQIL", "3": "MHKFRYW"}
# '1'stand for (0-2.78); '2'stand for (2.95-4.0), '3' stand for (4.03-8.08)

_Polarity = {"1": "LIFWCMVY", "2": "CPNVEQIL", "3": "KMHFRYW"}
# '1'stand for (4.9-6.2); '2'stand for (8.0-9.2), '3' stand for (10.4-13.0)

_Charge = {"1": "KR", "2": "ANCQGHILMFPSTWYV", "3": "DE"}
# '1'stand for Positive; '2'stand for Neutral, '3' stand for Negative

_SecondaryStr = {"1": "EALMQKRH", "2": "VIYCWFT", "3": "GNPSD"}
# '1'stand for Helix; '2'stand for Strand, '3' stand for coil

_SolventAccessibility = {"1": "ALFCGIVW", "2": "RKQEND", "3": "MPSTHY"}
# '1'stand for Buried; '2'stand for Exposed, '3' stand for Intermediate

_Polarizability = {"1": "GASDT", "2": "CPNVEQIL", "3": "KMHFRYW"}
# '1'stand for (0-0.108); '2'stand for (0.128-0.186), '3' stand for (0.219-0.409)


# You can continuely add other properties of AADs to compute descriptors of
# protein sequence.

_AATProperty = (
    _Hydrophobicity,
    _NormalizedVDWV,
    _Polarity,
    _Charge,
    _SecondaryStr,
    _SolventAccessibility,
    _Polarizability,
)

_AATPropertyName = (
    "_Hydrophobicity",
    "_NormalizedVDWV",
    "_Polarity",
    "_Charge",
    "_SecondaryStr",
    "_SolventAccessibility",
    "_Polarizability",
)


def StringtoNum(ProteinSequence, AAProperty):
    """
    Tranform the protein sequence into the string form such as 32123223132121123.

    Examples
    --------
    >>> result = StringtoNum(protein,AAProperty)

    Parameters
    ----------
     ProteinSequence : str
        a pure protein sequence

    AAProperty is a dict form containing classifciation of amino acids such as _Polarizability.

    Returns
    -------
     result is a string such as 123321222132111123222
    """
    hardProteinSequence = copy.deepcopy(ProteinSequence)
    for k, m in list(AAProperty.items()):
        for index in m:
            hardProteinSequence = hardProteinSequence.replace(index, k)
    TProteinSequence = hardProteinSequence

    return TProteinSequence


def CalculateComposition(ProteinSequence, AAProperty, AAPName):
    """
    A method used for computing composition descriptors.

    Examples
    --------
    >>> result = CalculateComposition(protein,AAProperty,AAPName)

    Parameters
    ----------
     ProteinSequence : str
        a pure protein sequence

    AAProperty is a dict form containing classifciation of amino acids such as _Polarizability.

    AAPName is a string used for indicating a AAP name.

    Returns
    -------
     result is a dict form containing composition descriptors based on the given property.
    """
    TProteinSequence = StringtoNum(ProteinSequence, AAProperty)
    Result = {}
    Num = len(TProteinSequence)
    Result[AAPName + "C" + "1"] = round(float(TProteinSequence.count("1")) / Num, 3)
    Result[AAPName + "C" + "2"] = round(float(TProteinSequence.count("2")) / Num, 3)
    Result[AAPName + "C" + "3"] = round(float(TProteinSequence.count("3")) / Num, 3)
    return Result


def CalculateTransition(ProteinSequence, AAProperty, AAPName):
    """
    A method used for computing transition descriptors

    Examples
    --------
    >>> result = CalculateTransition(protein,AAProperty,AAPName)

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    AAProperty is a dict form containing classifciation of amino acids such as _Polarizability.

    AAPName is a string used for indicating a AAP name.

    Returns
    -------
    result is a dict form containing transition descriptors based on the given property.
    """
    TProteinSequence = StringtoNum(ProteinSequence, AAProperty)
    Result = {}
    Num = len(TProteinSequence)
    CTD = TProteinSequence
    Result[AAPName + "T" + "12"] = round(
        float(CTD.count("12") + CTD.count("21")) / (Num - 1), 3
    )
    Result[AAPName + "T" + "13"] = round(
        float(CTD.count("13") + CTD.count("31")) / (Num - 1), 3
    )
    Result[AAPName + "T" + "23"] = round(
        float(CTD.count("23") + CTD.count("32")) / (Num - 1), 3
    )
    return Result


def CalculateDistribution(
    ProteinSequence: str, AAProperty: Dict[Any, Any], AAPName: str
):
    """
    A method used for computing distribution descriptors.

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence.
    AAProperty : Dict[Any, Any]
        contains classifciation of amino acids such as _Polarizability
    AAPName : str

    Returns
    -------
    result is a dict form containing Distribution descriptors based on the given property.

    Examples
    --------
    >>> result = CalculateDistribution(protein,AAProperty,AAPName)
    """
    TProteinSequence = StringtoNum(ProteinSequence, AAProperty)
    Result: Dict[str, float] = {}
    Num = len(TProteinSequence)
    for i in ("1", "2", "3"):
        num = TProteinSequence.count(i)
        ink = 1
        indexk = 0
        cds = []
        while ink <= num:
            indexk = TProteinSequence.find(i, indexk) + 1
            cds.append(indexk)
            ink = ink + 1

        if cds == []:
            Result[AAPName + "D" + i + "001"] = 0
            Result[AAPName + "D" + i + "025"] = 0
            Result[AAPName + "D" + i + "050"] = 0
            Result[AAPName + "D" + i + "075"] = 0
            Result[AAPName + "D" + i + "100"] = 0
        else:
            Result[AAPName + "D" + i + "001"] = round(float(cds[0]) / Num * 100, 3)
            Result[AAPName + "D" + i + "025"] = round(
                float(cds[int(math.floor(num * 0.25)) - 1]) / Num * 100, 3
            )
            Result[AAPName + "D" + i + "050"] = round(
                float(cds[int(math.floor(num * 0.5)) - 1]) / Num * 100, 3
            )
            Result[AAPName + "D" + i + "075"] = round(
                float(cds[int(math.floor(num * 0.75)) - 1]) / Num * 100, 3
            )
            Result[AAPName + "D" + i + "100"] = round(float(cds[-1]) / Num * 100, 3)

    return Result


def CalculateCompositionHydrophobicity(ProteinSequence: str):
    """
    A method used for calculating composition descriptors based on
    Hydrophobicity of AADs.

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    Returns
    -------
    result is a dict form containing Composition descriptors based on Hydrophobicity.

    Examples
    --------
    >>> result = CalculateCompositionHydrophobicity(protein)
    """
    result = CalculateComposition(ProteinSequence, _Hydrophobicity, "_Hydrophobicity")
    return result


def CalculateCompositionNormalizedVDWV(ProteinSequence: str):
    """
    A method used for calculating composition descriptors based on
    NormalizedVDWV of AADs.

    Examples
    --------
    >>> result = CalculateCompositionNormalizedVDWV(protein)

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    Returns
    -------
    result is a dict form containing Composition descriptors based on NormalizedVDWV.
    """
    result = CalculateComposition(ProteinSequence, _NormalizedVDWV, "_NormalizedVDWV")
    return result


def CalculateCompositionPolarity(ProteinSequence: str):
    """
    A method used for calculating composition descriptors based on Polarity of
    AADs.

    Examples
    --------
    >>> result = CalculateCompositionPolarity(protein)

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    Returns
    -------
    result is a dict form containing Composition descriptors based on Polarity.
    """
    result = CalculateComposition(ProteinSequence, _Polarity, "_Polarity")
    return result


def CalculateCompositionCharge(ProteinSequence: str):
    """
    A method used for calculating composition descriptors based on Charge of
    AADs.

    Examples
    --------
    >>> result = CalculateCompositionCharge(protein)

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    Returns
    -------
    result is a dict form containing Composition descriptors based on Charge.
    """
    result = CalculateComposition(ProteinSequence, _Charge, "_Charge")
    return result


def CalculateCompositionSecondaryStr(ProteinSequence: str):
    """
    A method used for calculating composition descriptors based on SecondaryStr
    of AADs.

    Examples
    --------
    >>> result = CalculateCompositionSecondaryStr(protein)

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    Returns
    -------
    result is a dict form containing Composition descriptors based on SecondaryStr.
    """
    result = CalculateComposition(ProteinSequence, _SecondaryStr, "_SecondaryStr")
    return result


def CalculateCompositionSolventAccessibility(ProteinSequence: str):
    """
    A method used for calculating composition descriptors based on SolventAccessibility
    of  AADs.

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    Returns
    -------
    result is a dict form containing Composition descriptors based on SolventAccessibility.

    Examples
    --------
    >>> result = CalculateCompositionSolventAccessibility(protein)
    """
    result = CalculateComposition(
        ProteinSequence, _SolventAccessibility, "_SolventAccessibility"
    )
    return result


def CalculateCompositionPolarizability(ProteinSequence: str):
    """
    A method used for calculating composition descriptors based on Polarizability of

    AADs.

    Examples
    --------
    >>> result = CalculateCompositionPolarizability(protein)

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    Returns
    -------
    result is a dict form containing Composition descriptors based on Polarizability.
    """
    result = CalculateComposition(ProteinSequence, _Polarizability, "_Polarizability")
    return result


def CalculateTransitionHydrophobicity(ProteinSequence: str):
    """
    A method used for calculating Transition descriptors based on Hydrophobicity of

    AADs.

    Examples
    --------
    >>> result = CalculateTransitionHydrophobicity(protein)

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    Returns
    -------
    result is a dict form containing Transition descriptors based on Hydrophobicity.
    """
    result = CalculateTransition(ProteinSequence, _Hydrophobicity, "_Hydrophobicity")
    return result


def CalculateTransitionNormalizedVDWV(ProteinSequence: str):
    """
    A method used for calculating Transition descriptors based on NormalizedVDWV of

    AADs.

    Examples
    --------
    >>> result = CalculateTransitionNormalizedVDWV(protein)

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    Returns
    -------
    result is a dict form containing Transition descriptors based on NormalizedVDWV.
    """
    result = CalculateTransition(ProteinSequence, _NormalizedVDWV, "_NormalizedVDWV")
    return result


def CalculateTransitionPolarity(ProteinSequence: str):
    """
    A method used for calculating Transition descriptors based on Polarity of
    AADs.

    Examples
    --------
    >>> result = CalculateTransitionPolarity(protein)

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    Returns
    -------
    result is a dict form containing Transition descriptors based on Polarity.
    """
    result = CalculateTransition(ProteinSequence, _Polarity, "_Polarity")
    return result


def CalculateTransitionCharge(ProteinSequence: str):
    """
    A method used for calculating Transition descriptors based on Charge of
    AADs.

    Examples
    --------
    >>> result = CalculateTransitionCharge(protein)

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    Returns
    -------
    result is a dict form containing Transition descriptors based on Charge.
    """
    result = CalculateTransition(ProteinSequence, _Charge, "_Charge")
    return result


def CalculateTransitionSecondaryStr(ProteinSequence: str):
    """
    A method used for calculating Transition descriptors based on SecondaryStr of

    AADs.

    Examples
    --------
    >>> result = CalculateTransitionSecondaryStr(protein)

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    Returns
    -------
    result is a dict form containing Transition descriptors based on SecondaryStr.
    """
    result = CalculateTransition(ProteinSequence, _SecondaryStr, "_SecondaryStr")
    return result


def CalculateTransitionSolventAccessibility(ProteinSequence: str):
    """
    A method used for calculating Transition descriptors based on SolventAccessibility

    of  AADs.

    Examples
    --------
    >>> result = CalculateTransitionSolventAccessibility(protein)

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    Returns
    -------
    result is a dict form containing Transition descriptors based on SolventAccessibility.
    """
    result = CalculateTransition(
        ProteinSequence, _SolventAccessibility, "_SolventAccessibility"
    )
    return result


def CalculateTransitionPolarizability(ProteinSequence: str):
    """
    A method used for calculating Transition descriptors based on Polarizability of

    AADs.

    Examples
    --------
    >>> result = CalculateTransitionPolarizability(protein)

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    Returns
    -------
    result is a dict form containing Transition descriptors based on Polarizability.
    """
    result = CalculateTransition(ProteinSequence, _Polarizability, "_Polarizability")
    return result


def CalculateDistributionHydrophobicity(ProteinSequence: str):
    """
    A method used for calculating Distribution descriptors based on Hydrophobicity of

    AADs.

    Examples
    --------
    >>> result = CalculateDistributionHydrophobicity(protein)

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    Returns
    -------
    result is a dict form containing Distribution descriptors based on Hydrophobicity.
    """
    result = CalculateDistribution(ProteinSequence, _Hydrophobicity, "_Hydrophobicity")
    return result


def CalculateDistributionNormalizedVDWV(ProteinSequence: str):
    """
    A method used for calculating Distribution descriptors based on NormalizedVDWV of

    AADs.

    Examples
    --------
    >>> result = CalculateDistributionNormalizedVDWV(protein)

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    Returns
    -------
    result is a dict form containing Distribution descriptors based on NormalizedVDWV.
    """
    result = CalculateDistribution(ProteinSequence, _NormalizedVDWV, "_NormalizedVDWV")
    return result


def CalculateDistributionPolarity(ProteinSequence: str):
    """
    A method used for calculating Distribution descriptors based on Polarity of
    AADs.

    Examples
    --------
    >>> result = CalculateDistributionPolarity(protein)

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    Returns
    -------
    result is a dict form containing Distribution descriptors based on Polarity.
    """
    result = CalculateDistribution(ProteinSequence, _Polarity, "_Polarity")
    return result


def CalculateDistributionCharge(ProteinSequence: str):
    """
    A method used for calculating Distribution descriptors based on Charge of
    AADs.

    Examples
    --------
    >>> result = CalculateDistributionCharge(protein)

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    Returns
    -------
    result is a dict form containing Distribution descriptors based on Charge.
    """
    result = CalculateDistribution(ProteinSequence, _Charge, "_Charge")
    return result


def CalculateDistributionSecondaryStr(ProteinSequence: str):
    """
    A method used for calculating Distribution descriptors based on
    SecondaryStr of AADs.

    Examples
    --------
    >>> result = CalculateDistributionSecondaryStr(protein)

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    Returns
    -------
    result is a dict form containing Distribution descriptors based on SecondaryStr.
    """
    result = CalculateDistribution(ProteinSequence, _SecondaryStr, "_SecondaryStr")
    return result


def CalculateDistributionSolventAccessibility(ProteinSequence: str):

    """
    A method used for calculating Distribution descriptors based on SolventAccessibility

    of  AADs.

    Examples
    --------
    >>> result = CalculateDistributionSolventAccessibility(protein)

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    Returns
    -------
    result is a dict form containing Distribution descriptors based on SolventAccessibility.
    """
    result = CalculateDistribution(
        ProteinSequence, _SolventAccessibility, "_SolventAccessibility"
    )
    return result


def CalculateDistributionPolarizability(ProteinSequence: str):
    """
    A method used for calculating Distribution descriptors based on Polarizability of

    AADs.

    Examples
    --------
    >>> result = CalculateDistributionPolarizability(protein)

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    Returns
    -------
    result is a dict form containing Distribution descriptors based on Polarizability.
    """
    result = CalculateDistribution(ProteinSequence, _Polarizability, "_Polarizability")
    return result


def CalculateC(ProteinSequence: str):
    """
    Calculate all composition descriptors based seven different properties of
    AADs.

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    Returns
    -------
    result is a dict form containing all composition descriptors.

    Examples
    --------
    >>> result = CalculateC(protein)
    """
    result: Dict[Any, Any] = {}
    result.update(CalculateCompositionPolarizability(ProteinSequence))
    result.update(CalculateCompositionSolventAccessibility(ProteinSequence))
    result.update(CalculateCompositionSecondaryStr(ProteinSequence))
    result.update(CalculateCompositionCharge(ProteinSequence))
    result.update(CalculateCompositionPolarity(ProteinSequence))
    result.update(CalculateCompositionNormalizedVDWV(ProteinSequence))
    result.update(CalculateCompositionHydrophobicity(ProteinSequence))
    return result


def CalculateT(ProteinSequence: str) -> Dict[Any, Any]:
    """
    Calculate all transition descriptors based seven different properties of AADs.

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    Returns
    -------
    result is a dict form containing all transition descriptors.

    Examples
    --------
    >>> result = CalculateT(protein)
    """
    result: Dict[Any, Any] = {}
    result.update(CalculateTransitionPolarizability(ProteinSequence))
    result.update(CalculateTransitionSolventAccessibility(ProteinSequence))
    result.update(CalculateTransitionSecondaryStr(ProteinSequence))
    result.update(CalculateTransitionCharge(ProteinSequence))
    result.update(CalculateTransitionPolarity(ProteinSequence))
    result.update(CalculateTransitionNormalizedVDWV(ProteinSequence))
    result.update(CalculateTransitionHydrophobicity(ProteinSequence))
    return result


def CalculateD(ProteinSequence: str):
    """
    Calculate all distribution descriptors based seven different properties of AADs.

    Examples
    --------
    >>> result = CalculateD(protein)

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    Returns
    -------
    result is a dict form containing all distribution descriptors.
    """
    result: Dict[Any, Any] = {}
    result.update(CalculateDistributionPolarizability(ProteinSequence))
    result.update(CalculateDistributionSolventAccessibility(ProteinSequence))
    result.update(CalculateDistributionSecondaryStr(ProteinSequence))
    result.update(CalculateDistributionCharge(ProteinSequence))
    result.update(CalculateDistributionPolarity(ProteinSequence))
    result.update(CalculateDistributionNormalizedVDWV(ProteinSequence))
    result.update(CalculateDistributionHydrophobicity(ProteinSequence))
    return result


def CalculateCTD(ProteinSequence: str):
    """
    Calculate all CTD descriptors based seven different properties of AADs.

    Examples
    --------
    >>> result = CalculateCTD(protein)

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    Returns
    -------
    result is a dict form containing all CTD descriptors.
    """
    result: Dict[Any, Any] = {}
    result.update(CalculateCompositionPolarizability(ProteinSequence))
    result.update(CalculateCompositionSolventAccessibility(ProteinSequence))
    result.update(CalculateCompositionSecondaryStr(ProteinSequence))
    result.update(CalculateCompositionCharge(ProteinSequence))
    result.update(CalculateCompositionPolarity(ProteinSequence))
    result.update(CalculateCompositionNormalizedVDWV(ProteinSequence))
    result.update(CalculateCompositionHydrophobicity(ProteinSequence))
    result.update(CalculateTransitionPolarizability(ProteinSequence))
    result.update(CalculateTransitionSolventAccessibility(ProteinSequence))
    result.update(CalculateTransitionSecondaryStr(ProteinSequence))
    result.update(CalculateTransitionCharge(ProteinSequence))
    result.update(CalculateTransitionPolarity(ProteinSequence))
    result.update(CalculateTransitionNormalizedVDWV(ProteinSequence))
    result.update(CalculateTransitionHydrophobicity(ProteinSequence))
    result.update(CalculateDistributionPolarizability(ProteinSequence))
    result.update(CalculateDistributionSolventAccessibility(ProteinSequence))
    result.update(CalculateDistributionSecondaryStr(ProteinSequence))
    result.update(CalculateDistributionCharge(ProteinSequence))
    result.update(CalculateDistributionPolarity(ProteinSequence))
    result.update(CalculateDistributionNormalizedVDWV(ProteinSequence))
    result.update(CalculateDistributionHydrophobicity(ProteinSequence))
    return result
