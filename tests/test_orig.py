# -*- coding: utf-8 -*-
"""
Test all commonly used functions of propy.

Authors: Dongsheng Cao and Yizeng Liang.
Date: 2012.09.04
Email: oriental-cds@163.com
"""


def test_original():
    # First party
    import propy.AAComposition as AAC
    import propy.Autocorrelation as AC
    import propy.CTD as CTD
    import propy.GetProteinFromUniprot as GPFU
    import propy.GetSubSeq as GSS
    import propy.PseudoAAC as PAAC
    import propy.QuasiSequenceOrder as QSO

    print("testing the GetProteinFromUniprot module")
    ProteinSequence = GPFU.GetProteinSequence("P08172")

    print("testing the GetSubSeq module")
    sub_sequence = GSS.GetSubSequence(ProteinSequence, ToAA="D", window=5)
    print(sub_sequence)

    print("testing the AAComposition module")
    aa_composition = AAC.CalculateAAComposition(ProteinSequence)
    print(aa_composition)

    # Just call it. Would be nice to know what the expected return value is
    AAC.CalculateDipeptideComposition(ProteinSequence)
    AAC.GetSpectrumDict(ProteinSequence)
    AAC.CalculateAADipeptideComposition(ProteinSequence)

    print("testing the Autocorrelation module")
    normalized_moreau_broto_auto = AC.CalculateNormalizedMoreauBrotoAuto(
        ProteinSequence, [AC._ResidueASA], ["ResidueASA"]
    )
    print(normalized_moreau_broto_auto)

    moran_auto = AC.CalculateMoranAuto(
        ProteinSequence, [AC._ResidueASA], ["ResidueASA"]
    )
    print(moran_auto)
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
