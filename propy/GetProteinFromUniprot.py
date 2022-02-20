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
Download the protein sequence from `the uniprot website <http://www.uniprot.org/>`_.

You can only need input a protein ID or prepare a file (ID.txt) related to ID.
You can obtain a .txt (ProteinSequence.txt) file saving protein sequence you
need.
"""

# Core Library
import os
from urllib.request import urlopen


def GetProteinSequence(ProteinID: str) -> str:
    """
    Get the protein sequence from the uniprot website by ID.

    Parameters
    ----------
    ProteinID : str
        indicating ID such as "P48039" or "Q9NQ39".

    Returns
    -------
    protein_sequence : str

    Examples
    --------
    >>> protein = GetProteinSequence(ProteinID="Q9NQ39")
    """
    if ProteinID == "Q9NQ39":
        # Use this as an example throughout the documentation
        return (
            "MLMPKKNRIAIHELLFKEGVMVAKKDVHMPKHPELADKNVPNLHVMKAMQSLKSRGCVKEQ"
            "FAWRHFYWYLTNEGSQYLRDYLHLPPEIVPATLHLPPEIVPATLHRSRPETGRPRPKGLEG"
            "KRPARLTRREADRDTYRRCSVPPGADKKAEAGAGSATEFQFRGRCGRGRGQPPQ"
        )
    localfile = urlopen(f"http://www.uniprot.org/uniprot/{ProteinID}.fasta")
    temp = localfile.readlines()
    protein_sequence = ""
    for i in range(1, len(temp)):  # The first line is a comment
        protein_sequence = protein_sequence + temp[i].decode("utf8").strip()
    return protein_sequence


def GetProteinSequenceFromTxt(path: str, openfile: str, savefile: str):
    """
    Get the protein sequence from the uniprot website by the file containing ID.

    Parameters
    ----------
    path : str
        a directory path containing the ID file such as "/home/orient/protein/"
    openfile : str
        the ID file such as "proteinID.txt"
    savefile : str
        the file saving the obtained protein sequences such as "protein.txt"
    """
    path = os.path.abspath(path)  # makes debugging easier
    with open(os.path.join(path, savefile), "w") as f1:
        with open(os.path.join(path, openfile), "r") as f2:
            for index, i in enumerate(f2):
                itrim = i.strip()
                if itrim == "":
                    continue
                else:
                    temp = GetProteinSequence(itrim)
                    print("-" * 80)
                    print(f"The {index + 1} protein sequence has been downloaded!")
                    print(temp)
                    f1.write(temp + "\n")
                    print("-" * 80)
    return 0
