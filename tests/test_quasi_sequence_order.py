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
