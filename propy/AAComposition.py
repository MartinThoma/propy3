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
The module is used for computing the composition of amino acids, dipetide and
3-mers (tri-peptide) for a given protein sequence.

References
----------
.. [1] Reczko, M. and Bohr, H. (1994) The DEF data base of sequence based protein
   fold class predictions. Nucleic Acids Res, 22, 3616-3619.

.. [2] Hua, S. and Sun, Z. (2001) Support vector machine approach for protein
   subcellular localization prediction. Bioinformatics, 17, 721-728.

.. [3] Grassmann, J., Reczko, M., Suhai, S. and Edler, L. (1999) Protein fold
   class prediction: new methods of statistical classification. Proc Int Conf
   Intell Syst Mol Biol, 106-112.
"""

# Core Library
import re
from typing import Any, Dict, List

# First party
from propy import AALetter


def CalculateAAComposition(ProteinSequence: str) -> Dict[str, float]:
    """
    Calculate the composition of Amino acids for a given protein sequence.

    Parameters
    ----------
    ProteinSequence: str
        a pure protein sequence

    Returns
    -------
    result : Dict[str, float]
        contains the composition of 20 amino acids.

    Examples
    --------
    >>> from propy.GetProteinFromUniprot import GetProteinSequence
    >>> protein = GetProteinSequence(ProteinID="Q9NQ39")
    >>> result = CalculateAAComposition(protein)
    """
    sequence_length = len(ProteinSequence)
    result: Dict[str, float] = {}
    for i in AALetter:
        result[i] = round(float(ProteinSequence.count(i)) / sequence_length * 100, 3)
    return result


def CalculateDipeptideComposition(ProteinSequence: str) -> Dict[str, float]:
    """
    Calculate the composition of dipeptidefor a given protein sequence.

    Parameters
    ----------
    ProteinSequence : a pure protein sequence

    Returns
    -------
    result : Dict[str, float]
        contains the composition of 400 dipeptides

    Examples
    --------
    >>> from propy.GetProteinFromUniprot import GetProteinSequence
    >>> protein = GetProteinSequence(ProteinID="Q9NQ39")
    >>> result = CalculateDipeptideComposition(protein)
    """
    sequence_length = len(ProteinSequence)
    result = {}
    for i in AALetter:
        for j in AALetter:
            dipeptide = i + j
            result[dipeptide] = round(
                float(ProteinSequence.count(dipeptide)) / (sequence_length - 1) * 100, 2
            )
    return result


def Getkmers() -> List[str]:
    """
    Get the amino acid list of 3-mers.

    Returns
    -------
    result : List[str]
        contains 8000 tri-peptides

    Examples
    --------
    >>> result = Getkmers()
    """
    kmers = []
    for i in AALetter:
        for j in AALetter:
            for k in AALetter:
                kmers.append(i + j + k)
    return kmers


def GetSpectrumDict(proteinsequence: str) -> Dict[str, int]:
    """
    Calcualte the spectrum descriptors of 3-mers for a given protein.

    Parameters
    ----------
    proteinsequence : a pure protein sequence

    Returns
    -------
    result : Dict[str, int]
        contains the composition values of 8000 3-mers

    Examples
    --------
    >>> from propy.GetProteinFromUniprot import GetProteinSequence
    >>> protein = GetProteinSequence(ProteinID="Q9NQ39")
    >>> result = GetSpectrumDict(protein)
    """
    result = {}
    kmers = Getkmers()
    for i in kmers:
        result[i] = len(re.findall(i, proteinsequence))
    return result


def CalculateAADipeptideComposition(ProteinSequence: str) -> Dict[str, float]:
    """
    Calculate the composition of AADs, dipeptide and 3-mers for a given protein
    sequence.

    Parameters
    ----------
    ProteinSequence : str
        a pure protein sequence

    Returns
    -------
    result : Dict[str, float]
        contains all composition values of AADs, dipeptide and 3-mers (8420).

    Examples
    --------
    >>> from propy.GetProteinFromUniprot import GetProteinSequence
    >>> protein = GetProteinSequence(ProteinID="Q9NQ39")
    >>> result = CalculateAADipeptideComposition(protein)
    """
    result: Dict[Any, Any] = {}
    result.update(CalculateAAComposition(ProteinSequence))
    result.update(CalculateDipeptideComposition(ProteinSequence))
    result.update(GetSpectrumDict(ProteinSequence))
    return result
