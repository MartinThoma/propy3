# Third party
import pytest

# First party
from propy import PyPro
from propy.GetProteinFromUniprot import GetProteinSequence as gps


def test_docs():
    uniprotid = "P48039"
    proseq = gps(uniprotid)


@pytest.mark.xfail()
def test_p33765():
    proteinsequence = gps("P33765")  # download the protein sequence by uniprot id
    DesObject = PyPro.GetProDes(proteinsequence)  # construct a GetProDes object
    print(DesObject.GetCTD())  # calculate 147 CTD descriptors
    print(DesObject.GetAAComp())  # calculate 20 amino acid composition descriptors
    paac = DesObject.GetPAAC(
        lamda=10, weight=0.05
    )  # calculate 30 pseudo amino acid composition descriptors

    for i in paac:
        print(i, paaci)
