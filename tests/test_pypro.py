# Third party
import pytest


@pytest.mark.skip(reason="Currently fails on travis-ci.org with timeout")
def test_main():
    from propy.PyPro import GetProDes
    from propy.Autocorrelation import _Steric
    from propy.PseudoAAC import _Hydrophobicity, _hydrophilicity

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
