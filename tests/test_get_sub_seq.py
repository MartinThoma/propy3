# First party
from propy.GetSubSeq import GetSubSequence


def test_main():
    protein = "ADGCGVGEGTGQGPMCNCMCMKWVYADEDAADLESDSFADEDASLESDSFPWSNQRVFCSFADEDAS"
    subseq = GetSubSequence(protein, ToAA="D", window=10)
    print(subseq)
    print(len(subseq))
    # print(len(subseq[0]))
