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
