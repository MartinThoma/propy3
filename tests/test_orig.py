# -*- coding: utf-8 -*-
"""
Test all commonly used functions of propy.

Authors: Dongsheng Cao and Yizeng Liang.
Date: 2012.09.04
Email: oriental-cds@163.com
"""


def test_original():
    import propy.AAComposition as AAC
    import propy.Autocorrelation as AC
    import propy.CTD as CTD
    import propy.QuasiSequenceOrder as QSO
    import propy.PseudoAAC as PAAC
    import propy.GetProteinFromUniprot as GPFU
    import propy.GetSubSeq as GSS

    print("testing the GetProteinFromUniprot module")
    ProteinSequence = GPFU.GetProteinSequence("P08172")

    print("testing the GetSubSeq module")
    temp = GSS.GetSubSequence(ProteinSequence, ToAA="D", window=5)
    print(temp)

    print("testing the AAComposition module")
    temp = AAC.CalculateAAComposition(ProteinSequence)
    print(temp)

    temp = AAC.CalculateDipeptideComposition(ProteinSequence)
    temp = AAC.GetSpectrumDict(ProteinSequence)
    temp = AAC.CalculateAADipeptideComposition(ProteinSequence)

    print("testing the Autocorrelation module")
    temp = AC.CalculateNormalizedMoreauBrotoAuto(
        ProteinSequence, [AC._ResidueASA], ["ResidueASA"]
    )
    print(temp)

    temp = AC.CalculateMoranAuto(ProteinSequence, [AC._ResidueASA], ["ResidueASA"])
    print(temp)
    temp = AC.CalculateGearyAuto(ProteinSequence, [AC._ResidueASA], ["ResidueASA"])
    print(temp)
    temp = AC.CalculateAutoTotal(ProteinSequence)

    print("testing the CTD module")
    temp = CTD.CalculateC(ProteinSequence)
    print(temp)
    temp = CTD.CalculateT(ProteinSequence)
    print(temp)
    temp = CTD.CalculateD(ProteinSequence)
    print(temp)
    temp = CTD.CalculateCTD(ProteinSequence)
    print(temp)

    print("...............................................................")
    print("testing the QuasiSequenceOrder module")
    temp = QSO.GetSequenceOrderCouplingNumberTotal(ProteinSequence, maxlag=30)
    print(temp)
    temp = QSO.GetQuasiSequenceOrder(ProteinSequence, maxlag=30, weight=0.1)
    print(temp)

    print("...............................................................")
    print("testing the PseudoAAC module")
    temp = PAAC.GetAPseudoAAC(ProteinSequence, lamda=10, weight=0.5)
    print(temp)
    temp = PAAC._GetPseudoAAC(ProteinSequence, lamda=10, weight=0.05)
    print(temp)

    print("...............................................................")
    print("Tested successfully!")
