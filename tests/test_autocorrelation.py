# First party
from propy.Autocorrelation import (
    CalculateAutoTotal,
    CalculateMoranAutoMutability,
    CalculateNormalizedMoreauBrotoAuto,
    _AAProperty,
    _AAPropertyName,
)


def test_main():
    protein = "ADGCGVGEGTGQGPMCNCMCMKWVYADEDAADLESDSFADEDASLESDSFPWSNQRVFCSFADEDAS"
    temp1 = CalculateNormalizedMoreauBrotoAuto(
        protein, AAProperty=_AAProperty, AAPropertyName=_AAPropertyName
    )
    print(temp1)
    temp2 = CalculateMoranAutoMutability(protein)
    print(temp2)
    print(len(CalculateAutoTotal(protein)))
