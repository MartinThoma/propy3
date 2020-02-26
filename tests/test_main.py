# First party
from propy.GetProteinFromUniprot import GetProteinSequence as gps


def test_docs():
    uniprotid = "P48039"
    proseq = gps(uniprotid)
