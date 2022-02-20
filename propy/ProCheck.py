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
Check whether the input protein sequence is a valid amino acid sequence.
"""

# First party
from propy import AALetter


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
    >>> from propy.GetProteinFromUniprot import GetProteinSequence
    >>> protein = GetProteinSequence(ProteinID="Q9NQ39")
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
