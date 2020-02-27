# First party
from propy.ProCheck import ProteinCheck


def test_main():
    protein = "ADGCGVGEGTGQGPMCNCMCMKWVYADEDAADLESDSFADEDASLESDSFPWSNQRVFCSFADEDASU"
    print(ProteinCheck(protein))
