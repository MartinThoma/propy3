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
# Third party
import pytest


@pytest.mark.skip(reason="Currently fails on travis-ci.org with timeout")
def test_main():
    # First party
    from propy.Autocorrelation import _Steric
    from propy.PseudoAAC import _hydrophilicity, _Hydrophobicity
    from propy.PyPro import GetProDes

    protein = "ADGCGVGEGTGQGPMCNCMCMKWVYADEDAADLESDSFADEDASLESDSFPWSNQRVFCSFADEDAS"
    cds = GetProDes(protein)

    # print cds.GetAAComp()
    # print cds.GetDPComp()
    # print cds.GetTPComp()
    # print cds.GetCTD()
    print(cds.GetPAAC(lamda=5))
    # print cds.GetALL()
    print(cds.GetMoreauBrotoAutop(AAP=_Steric, AAPName="Steric"))
    print(cds.GetMoranAutop(AAP=_Steric, AAPName="Steric"))
    print(cds.GetGearyAutop(AAP=_Steric, AAPName="Steric"))
    print(cds.GetPAACp(lamda=5, weight=0.05, AAP=[_Hydrophobicity, _hydrophilicity]))
    print(cds.GetSubSeq(ToAA="D", window=5))

    proper = cds.GetAAindex23("GRAR740104", path=None)
    # print cds.GetAAindex1('KRIW790103',path='/home/orient')

    print(cds.GetQSOp(maxlag=30, weight=0.1, distancematrix=proper))
    print(cds.GetSOCNp(maxlag=30, distancematrix=proper))
