# -*- coding: utf-8 -*-
"""
Download the protein sequence from [the uniprot website](http://www.uniprot.org/).

You can only need input a protein ID or prepare a file (ID.txt) related to ID.
You can obtain a .txt (ProteinSequence.txt) file saving protein sequence you
need.

Authors: Dongsheng Cao and Yizeng Liang.
Date: 2012.9.3
Email: oriental-cds@163.com
"""

# Core Library
import os
import urllib.error
import urllib.parse
import urllib.request


def GetProteinSequence(ProteinID: str):
    """
    Get the protein sequence from the uniprot website by ID.

    Parameters
    ----------
    ProteinID : str
        indicating ID such as "P48039".

    Returns
    -------
    result :
        a protein sequence.

    Examples
    --------
    >>> result = GetProteinSequence(ProteinID)
    """

    ID = str(ProteinID)
    localfile = urllib.request.urlopen(
        "http://www.uniprot.org/uniprot/" + ID + ".fasta"
    )
    temp = localfile.readlines()
    res = ""
    for i in range(1, len(temp)):
        res = res + temp[i].decode("utf8").strip()
    return res


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

    Examples
    --------
    >>> result = GetProteinSequenceFromTxt(path, openfile, savefile)
    """
    with open(os.path.join(path, savefile), "wb") as f1:
        with open(os.path.join(path, openfile), "r") as f2:
            for index, i in enumerate(f2):
                itrim = i.strip()
                if itrim == "":
                    continue
                else:
                    temp = GetProteinSequence(itrim)
                    print("-" * 80)
                    print("The {index + 1} protein sequence has been downloaded!")
                    print(temp)
                    f1.write(temp + "\n")
                    print("-" * 80)
    return 0
