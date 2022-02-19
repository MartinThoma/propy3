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
from propy.CTD import CalculateCTD


def test_main():
    # import scipy,string

    # result=scipy.zeros((268,147))
    # f=file('protein1.txt','r')
    # for i,j in enumerate(f):
    #     temp=CalculateCTD(j.strip())
    #     result[i,:]=temp.values()
    # scipy.savetxt('ResultNCTRER.txt', result, fmt='%15.5f',delimiter='')
    protein = "ADGCGVGEGTGQGPMCNCMCMKWVYADEDAADLESDSFADEDASLESDSFPWSNQRVFCSFADEDAS"
    # print StringtoNum(protein,_Hydrophobicity)
    # print CalculateComposition(protein,_Hydrophobicity,'_Hydrophobicity')
    # print CalculateTransition(protein,_Hydrophobicity,'_Hydrophobicity')
    # print CalculateDistribution(protein,_Hydrophobicity,'_Hydrophobicity')
    # print CalculateDistributionSolventAccessibility(protein)
    # print len(CalculateCTD(protein))
    # print len(CalculateC(protein))
    # print len(CalculateT(protein))
    # print len(CalculateD(protein))
    print(CalculateCTD(protein))
