# propy3, formerly protpy, is a Python package to compute protein descriptors
# Copyright (C) 2012 Dongsheng Cao and Yizeng Liang, oriental-cds@163.com
# Copyright (C) 2020-2022 Martin Thoma

# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor,
# Boston, MA  02110-1301, USA.
# Core Library
import os
from tempfile import mkstemp

# First party
from propy.GetProteinFromUniprot import GetProteinSequence, GetProteinSequenceFromTxt


def test_main():
    _, result_filepath = mkstemp(suffix="result.txt", prefix="propy3")
    _, target_filepath = mkstemp(suffix="target.txt", prefix="propy3")
    with open(result_filepath, "wb") as savefile:
        with open(target_filepath, "r") as localfile:
            for index, i in enumerate(localfile):
                itrim = i.strip()
                if itrim == "":
                    continue
                else:
                    temp = GetProteinSequence(itrim)
                    print("--------------------------------------------------------")
                    print("The %d protein sequence has been downloaded!" % (index + 1))
                    print(temp)
                    savefile.write((temp + "\n").encode("utf8"))
                    print("--------------------------------------------------------")

    flag = GetProteinSequenceFromTxt(
        "/home/orient/ProPy/", target_filepath, result_filepath
    )
    print(flag)

    # Cleanup
    os.remove(result_filepath)
    os.remove(target_filepath)
