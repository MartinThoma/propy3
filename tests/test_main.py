# propy3, formerly protpy, is a Python package to compute protein descriptors
# Copyright (C) 2012 Dongsheng Cao and Yizeng Liang, oriental-cds@163.com
# Copyright (C) 2020-2022 Martin Thoma

# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; in version 2
# of the License.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor,
# Boston, MA  02110-1301, USA.
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
def test_p33765():
    # TODO: "P33765" gives "HTTP Error 300" (abiguity?) Why?
    proteinsequence = gps("P33765")  # download the protein sequence by uniprot id
    DesObject = PyPro.GetProDes(proteinsequence)  # construct a GetProDes object
    print(DesObject.GetCTD())  # calculate 147 CTD descriptors
    print(DesObject.GetAAComp())  # calculate 20 amino acid composition descriptors
    paac = DesObject.GetPAAC(
        lamda=10, weight=0.05
    )  # calculate 30 pseudo amino acid composition descriptors

    for i in paac:
        print(i)
