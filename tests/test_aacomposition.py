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
# First party
from propy.AAComposition import (
    CalculateAAComposition,
    CalculateAADipeptideComposition,
    CalculateDipeptideComposition,
    GetSpectrumDict,
)


def test_main():
    protein = "ADGCGVGEGTGQGPMCNCMCMKWVYADEDAADLESDSFADEDASLESDSFPWSNQRVFCSFADEDAS"

    AAC = CalculateAAComposition(protein)
    print(AAC)
    DIP = CalculateDipeptideComposition(protein)
    print(DIP)
    spectrum = GetSpectrumDict(protein)
    print(spectrum)
    res = CalculateAADipeptideComposition(protein)
    print(len(res))
