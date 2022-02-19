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
from propy.QuasiSequenceOrder import GetSequenceOrderCouplingNumberTotal


def test_main():
    # from propy.QuasiSequenceOrder import GetQuasiSequenceOrderp

    protein = "ELRLRYCAPAGFALLKCNDADYDGFKTNCSNVSVVHCTNLMNTTVTTGLLLNGSYSENRT\
QIWQKHRTSNDSALILLNKHYNLTVTCKRPGNKTVLPVTIMAGLVFHSQKYNLRLRQAWC\
HFPSNWKGAWKEVKEEIVNLPKERYRGTNDPKRIFFQRQWGDPETANLWFNCHGEFFYCK\
MDWFLNYLNNLTVDADHNECKNTSGTKSGNKRAPGPCVQRTYVACHIRSVIIWLETISKK\
TYAPPREGHLECTSTVTGMTVELNYIPKNRTNVTLSPQIESIWAAELDRYKLVEITPIGF\
APTEVRRYTGGHERQKRVPFVVQSQHLLAGILQQQKNLLAAVEAQQQMLKLTIWGVK"
    print(len(protein))
    SCN = GetSequenceOrderCouplingNumberTotal(protein, maxlag=30)
    print(len(SCN))
    for i in SCN:
        print(i, SCN[i])
    #
    #     QSO1=GetQuasiSequenceOrder1(protein,maxlag=30,weight=0.1)
    #     print QSO1
    #     for i in QSO1:
    #         print i, QSO1[i]
    #
    #     QSO2=GetQuasiSequenceOrder2(protein,maxlag=30,weight=0.1)
    #     print QSO2
    #     for i in QSO2:
    #         print i, QSO2[i]
    #     QSO=GetQuasiSequenceOrder(protein,maxlag=30,weight=0.1)
    #     print len(QSO)
    #     for i in QSO:
    #         print i, QSO[i]

    #     SCN=GetSequenceOrderCouplingNumberp(protein,maxlag=30,distancematrix=_Distance1)
    #     print len(SCN)
    #     for i in SCN:
    #         print i, SCN[i]

    # QSO = GetQuasiSequenceOrderp(protein, maxlag=30, distancematrix=_Distance1)
    # print(len(QSO))
    # for i in QSO:
    #     print(i, QSO[i])
