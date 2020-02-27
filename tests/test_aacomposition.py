# First party
from propy.AAComposition import (
    CalculateAAComposition,
    CalculateAADipeptideComposition,
    CalculateDipeptideComposition,
    GetSpectrumDict,
)


def test_main():
    protein = "ADGCGVGEGTGQGPMCNCMCMKWVYADEDAADLESDSFADEDASLESDSFPWSNQRVFCSFADEDAS"

    AAC = CalculateAAComposition(protein)
    print(AAC)
    DIP = CalculateDipeptideComposition(protein)
    print(DIP)
    spectrum = GetSpectrumDict(protein)
    print(spectrum)
    res = CalculateAADipeptideComposition(protein)
    print(len(res))
