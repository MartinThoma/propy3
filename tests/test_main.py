# Third party
import pytest

# First party
from propy import PyPro
from propy.GetProteinFromUniprot import GetProteinSequence as gps


def test_docs():
    uniprotid = "P48039"
    gps(uniprotid)  # Check the return value!


def test_marina():
    # First party
    from propy import CTD
    from propy import AAComposition as AAC
    from propy.PyPro import GetProDes

    # Protein sequence
    proseq = (
        "MENATLLKSTTRHIRIFAAEIDRDGELVPSNQVLTLDIDPDNEFNWNEDALQKIYRKFDELV"
        "EASSGADLTDYNLRRIGSDLEHYLRSLLQKGEISYNLSARVTNYSLGLPQVAVEDK"
    )
    _ = AAC.CalculateAAComposition(proseq)  # TODO: Check return value

    _ = CTD.CalculateC(proseq)  # TODO: Check return value

    Des = GetProDes(proseq)
    alldes = Des.GetALL()
    for desc in alldes:
        print(desc, alldes[desc])


@pytest.mark.xfail()
def test_p48039():
    proteinsequence = gps("P48039")  # download the protein sequence by uniprot id
    DesObject = PyPro.GetProDes(proteinsequence)  # construct a GetProDes object
    print(DesObject.GetCTD())  # calculate 147 CTD descriptors
    print(DesObject.GetAAComp())  # calculate 20 amino acid composition descriptors
    paac = DesObject.GetPAAC(
        lamda=10, weight=0.05
    )  # calculate 30 pseudo amino acid composition descriptors

    for i in paac:
        print(i)
