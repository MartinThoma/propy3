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
