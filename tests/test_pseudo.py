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
from propy.PseudoAAC import GetPseudoAAC, _hydrophilicity, _Hydrophobicity


def test_main():
    protein = "MTDRARLRLHDTAAGVVRDFVPLRPGHVSIYLCGATVQGLPHIGHVRSGVAFDILRRWLL\
ARGYDVAFIRNVTDIEDKILAKAAAAGRPWWEWAATHERAFTAAYDALDVLPPSAEPRAT\
GHITQMIEMIERLIQAGHAYTGGGDVYFDVLSYPEYGQLSGHKIDDVHQGEGVAAGKRDQ\
RDFTLWKGEKPGEPSWPTPWGRGRPGWHLECSAMARSYLGPEFDIHCGGMDLVFPHHENE\
IAQSRAAGDGFARYWLHNGWVTMGGEKMSKSLGNVLSMPAMLQRVRPAELRYYLGSAHYR\
SMLEFSETAMQDAVKAYVGLEDFLHRVRTRVGAVCPGDPTPRFAEALDDDLSVPIALAEI\
HHVRAEGNRALDAGDHDGALRSASAIRAMMGILGCDPLDQRWESRDETSAALAAVDVLVQ\
AELQNREKAREQRNWALADEIRGRLKRAGIEVTDTADGPQWSLLGGDTK"
    protein = protein.strip()
    #     temp=_GetCorrelationFunction('S','D')
    #     print temp
    #
    #     print _GetSequenceOrderCorrelationFactor(protein,k=4)
    #
    #     PAAC1=_GetPseudoAAC1(protein,lamda=4)
    #     for i in PAAC1:
    #         print i, PAAC1[i]
    #     PAAC2=_GetPseudoAAC2(protein,lamda=4)
    #     for i in PAAC2:
    #         print i, PAAC2[i]
    #     print len(PAAC1)
    #     print _GetSequenceOrderCorrelationFactorForAPAAC(protein,k=1)
    #     APAAC1=_GetAPseudoAAC1(protein,lamda=4)
    #     for i in APAAC1:
    #         print i, APAAC1[i]

    #     APAAC2=GetAPseudoAAC2(protein,lamda=4)
    #     for i in APAAC2:
    #         print i, APAAC2[i]
    #     APAAC=GetAPseudoAAC(protein,lamda=4)
    #
    #     for i in APAAC:
    #         print i, APAAC[i]

    PAAC = GetPseudoAAC(protein, lamda=5, AAP=[_Hydrophobicity, _hydrophilicity])

    for i in PAAC:
        print(i, PAAC[i])
