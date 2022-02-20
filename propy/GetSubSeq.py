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
The prediction of functional sites (e.g. methylation) of proteins usually needs
to split the total protein into a set of segments around specific amino acid.
Given a specific window size p, we can obtain all segments of length equal to
(2*p+1) very easily. Note that the output of the method is a list form.
"""

# Core Library
import re
from typing import List

# First party
from propy import AALetter


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
    >>> from propy.GetProteinFromUniprot import GetProteinSequence
    >>> protein = GetProteinSequence(ProteinID="Q9NQ39")
    >>> result = GetSubSequence(protein)
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
