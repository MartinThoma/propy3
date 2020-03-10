# -*- coding: utf-8 -*-
"""
Download the protein sequence from `the uniprot website <http://www.uniprot.org/>`_.

You can only need input a protein ID or prepare a file (ID.txt) related to ID.
You can obtain a .txt (ProteinSequence.txt) file saving protein sequence you
need.

Authors: Dongsheng Cao and Yizeng Liang.
Date: 2012.9.3
Email: oriental-cds@163.com
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
