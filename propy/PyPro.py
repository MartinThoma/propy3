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
"""Computing different types of protein descriptors."""

# Core Library
from typing import Any, Dict, List, Optional

# Local
from .AAComposition import (
    CalculateAAComposition,
    CalculateDipeptideComposition,
    GetSpectrumDict,
)
from .AAIndex import GetAAIndex1, GetAAIndex23
from .Autocorrelation import (
    CalculateEachGearyAuto,
    CalculateEachMoranAuto,
    CalculateEachNormalizedMoreauBrotoAuto,
    CalculateGearyAutoTotal,
    CalculateMoranAutoTotal,
    CalculateNormalizedMoreauBrotoAutoTotal,
)
from .CTD import CalculateCTD
from .GetSubSeq import GetSubSequence
from .PseudoAAC import GetAPseudoAAC, GetPseudoAAC, _GetPseudoAAC
from .QuasiSequenceOrder import (
    GetQuasiSequenceOrder,
    GetQuasiSequenceOrderp,
    GetSequenceOrderCouplingNumberp,
    GetSequenceOrderCouplingNumberTotal,
)


class GetProDes:
    """Collect all descriptor calcualtion modules."""

    AALetter = list("ARNDCEQGHILKMFPSTWYV")

    Version = 1.0

    def __init__(self, ProteinSequence: str = "") -> None:
        """Input a protein sequence."""
        if len(ProteinSequence) == 0:
            print(
                "You must input a protein sequence "
                "when constructing a object. It is a string!"
            )
        else:
            self.ProteinSequence = ProteinSequence

    def GetAAComp(self) -> Dict[str, float]:
        """
        Amino acid compositon descriptors (20).

        Examples
        --------
        >>> from propy.GetProteinFromUniprot import GetProteinSequence
        >>> protein = GetProteinSequence(ProteinID="Q9NQ39")
        >>> result = GetProDes(protein).GetAAComp()
        """
        res = CalculateAAComposition(self.ProteinSequence)
        return res

    def GetDPComp(self) -> Dict[str, float]:
        """
        Dipeptide composition descriptors (400).

        Examples
        --------
        >>> from propy.GetProteinFromUniprot import GetProteinSequence
        >>> protein = GetProteinSequence(ProteinID="Q9NQ39")
        >>> result = GetProDes(protein).GetDPComp()
        """
        res = CalculateDipeptideComposition(self.ProteinSequence)
        return res

    def GetTPComp(self) -> Dict[str, int]:
        """
        Tri-peptide composition descriptors (8000).

        Examples
        --------
        >>> from propy.GetProteinFromUniprot import GetProteinSequence
        >>> protein = GetProteinSequence(ProteinID="Q9NQ39")
        >>> result = GetProDes(protein).GetTPComp()
        """
        res = GetSpectrumDict(self.ProteinSequence)
        return res

    def GetMoreauBrotoAuto(self) -> Dict[Any, Any]:
        """
        Normalized Moreau-Broto autocorrelation descriptors (240).

        Examples
        --------
        >>> from propy.GetProteinFromUniprot import GetProteinSequence
        >>> protein = GetProteinSequence(ProteinID="Q9NQ39")
        >>> result = GetProDes(protein).GetMoreauBrotoAuto()
        """
        res = CalculateNormalizedMoreauBrotoAutoTotal(self.ProteinSequence)
        return res

    def GetMoranAuto(self) -> Dict[Any, Any]:
        """
        Moran autocorrelation descriptors (240).

        Examples
        --------
        >>> from propy.GetProteinFromUniprot import GetProteinSequence
        >>> protein = GetProteinSequence(ProteinID="Q9NQ39")
        >>> result = GetProDes(protein).GetMoranAuto()
        """
        res = CalculateMoranAutoTotal(self.ProteinSequence)
        return res

    def GetGearyAuto(self) -> Dict[Any, Any]:
        """
        Geary autocorrelation descriptors (240).

        Examples
        --------
        >>> from propy.GetProteinFromUniprot import GetProteinSequence
        >>> protein = GetProteinSequence(ProteinID="Q9NQ39")
        >>> result = GetProDes(protein).GetGearyAuto()
        """
        res = CalculateGearyAutoTotal(self.ProteinSequence)
        return res

    def GetCTD(self) -> Dict[Any, Any]:
        """
        Composition Transition Distribution descriptors (147).

        Examples
        --------
        >>> from propy.GetProteinFromUniprot import GetProteinSequence
        >>> protein = GetProteinSequence(ProteinID="Q9NQ39")
        >>> result = GetProDes(protein).GetCTD()
        """
        res = CalculateCTD(self.ProteinSequence)
        return res

    def GetPAAC(self, lamda: int = 10, weight: float = 0.05) -> Dict[Any, Any]:
        """
        Type I Pseudo amino acid composition descriptors (default is 30).

        Parameters
        ----------
        lamda : int, optional (default: 10)
            reflects the rank of correlation and is a non-Negative integer,
            such as 15. Note that (1)lamda should NOT be larger than the length
            of input protein sequence; (2) lamda must be non-Negative integer,
            such as 0, 1, 2, ...; (3) when lamda =0, the output of PseAA server
            is the 20-D amino acid composition.
        weight : float, optional (default: 0.05)
            is designed for the users to put weight on the additional PseAA
            components with respect to the conventional AA components. The user
            can select any value within the region from 0.05 to 0.7 for the
            weight factor.

        Examples
        --------
        >>> from propy.GetProteinFromUniprot import GetProteinSequence
        >>> protein = GetProteinSequence(ProteinID="Q9NQ39")
        >>> result = GetProDes(protein).GetPAAC(lamda=10, weight=0.05)
        """
        res = _GetPseudoAAC(self.ProteinSequence, lamda=lamda, weight=weight)
        return res

    def GetPAACp(
        self, lamda: int = 10, weight: float = 0.05, AAP: Optional[List[Any]] = None
    ) -> Dict[Any, Any]:
        """
        Type I Pseudo amino acid composition descriptors for the given properties

        Default is 30.

        Parameters
        ----------
        lamda : int, optional (default: 10)
            reflects the rank of correlation and is a non-Negative integer,
            such as 15. Note that (1)lamda should NOT be larger than the length
            of input protein sequence; (2) lamda must be non-Negative integer,
            such as 0, 1, 2, ...; (3) when lamda =0, the output of PseAA server
            is the 20-D amino acid composition.
        weight : float, optional (default: 0.05)
            is designed for the users to put weight on the additional PseAA
            components with respect to the conventional AA components. The user
            can select any value within the region from 0.05 to 0.7 for the
            weight factor.
        AAP : List
            contains the properties, each of which is a dict form.

        Examples
        --------
        >>> from propy.GetProteinFromUniprot import GetProteinSequence
        >>> protein = GetProteinSequence(ProteinID="Q9NQ39")
        >>> result = GetProDes(protein).GetPAACp(lamda=10, weight=0.05, AAP=[])
        """
        if AAP is None:
            AAP = []
        res = GetPseudoAAC(self.ProteinSequence, lamda=lamda, weight=weight, AAP=AAP)
        return res

    def GetAPAAC(self, lamda: int = 10, weight: float = 0.5) -> Dict[Any, Any]:
        """
        Amphiphilic (Type II) Pseudo amino acid composition descriptors.

        default is 30

        Parameters
        ----------
        lamda : int, optional (default: 10)
            reflects the rank of correlation and is a non-Negative integer,
            such as 15. Note that (1)lamda should NOT be larger than the length
            of input protein sequence; (2) lamda must be non-Negative integer,
            such as 0, 1, 2, ...; (3) when lamda =0, the output of PseAA server
            is the 20-D amino acid composition.
        weight : float, optional (default: 0.05)
            is designed for the users to put weight on the additional PseAA
            components with respect to the conventional AA components. The user
            can select any value within the region from 0.05 to 0.7 for the
            weight factor.

        Examples
        --------
        >>> from propy.GetProteinFromUniprot import GetProteinSequence
        >>> protein = GetProteinSequence(ProteinID="Q9NQ39")
        >>> result = GetProDes(protein).GetAPAAC(lamda=10, weight=0.5)
        """
        res = GetAPseudoAAC(self.ProteinSequence, lamda=lamda, weight=weight)
        return res

    def GetSOCN(self, maxlag: int = 45) -> Dict[Any, Any]:
        """
        Sequence order coupling numbers  default is 45.

        Parameters
        ----------
        maxlag : int
            is the maximum lag and the length of the protein should be larger
            than maxlag

        Examples
        --------
        >>> from propy.GetProteinFromUniprot import GetProteinSequence
        >>> protein = GetProteinSequence(ProteinID="Q9NQ39")
        >>> result = GetProDes(protein).GetSOCN(maxlag=45)
        """
        res = GetSequenceOrderCouplingNumberTotal(self.ProteinSequence, maxlag=maxlag)
        return res

    def GetSOCNp(
        self, maxlag: int = 45, distancematrix: Optional[Dict[Any, Any]] = None
    ) -> Dict[Any, Any]:
        """
        Sequence order coupling numbers  default is 45.

        Parameters
        ----------
        maxlag is the maximum lag and the length of the protein should be larger
        than maxlag. default is 45.

        distancematrix is a dict form containing 400 distance values

        Examples
        --------
        >>> from propy.GetProteinFromUniprot import GetProteinSequence
        >>> protein = GetProteinSequence(ProteinID="Q9NQ39")
        >>> result = GetProDes(protein).GetSOCN(maxlag=45)
        """
        if distancematrix is None:
            distancematrix = {}
        res = GetSequenceOrderCouplingNumberp(
            self.ProteinSequence, maxlag=maxlag, distancematrix=distancematrix
        )
        return res

    def GetQSO(self, maxlag: int = 30, weight: float = 0.1) -> Dict[Any, Any]:
        """
        Quasi sequence order descriptors  default is 50.

        Parameters
        ----------
        result = GetQSO(maxlag=30, weight=0.1)

        maxlag is the maximum lag and the length of the protein should be larger

        than maxlag. default is 45.
        """
        res = GetQuasiSequenceOrder(self.ProteinSequence, maxlag=maxlag, weight=weight)
        return res

    def GetQSOp(
        self,
        maxlag: int = 30,
        weight: float = 0.1,
        distancematrix: Optional[Dict[Any, Any]] = None,
    ) -> Dict[Any, Any]:
        """
        Quasi sequence order descriptors  default is 50.

        Parameters
        ----------
        result = GetQSO(maxlag=30, weight=0.1)

        maxlag is the maximum lag and the length of the protein should be larger
        than maxlag. default is 45.

        distancematrix is a dict form containing 400 distance values
        """
        if distancematrix is None:
            distancematrix = {}
        res = GetQuasiSequenceOrderp(
            self.ProteinSequence,
            maxlag=maxlag,
            weight=weight,
            distancematrix=distancematrix,
        )
        return res

    def GetMoreauBrotoAutop(
        self, AAP: Optional[Dict[Any, Any]] = None, AAPName: str = "p"
    ) -> Dict[str, float]:
        """
        Normalized Moreau-Broto autocorrelation descriptors for the given property (30).

        Parameters
        ----------
        AAP : Dict[Any, Any]
            contains physicochemical properities of 20 amino acids

        Examples
        --------
        >>> from propy.GetProteinFromUniprot import GetProteinSequence
        >>> protein = GetProteinSequence(ProteinID="Q9NQ39")
        >>> result = GetProDes(protein).GetMoreauBrotoAutop(AAP={}, AAPName='p')
        """
        if AAP is None:
            AAP = {}
        res = CalculateEachNormalizedMoreauBrotoAuto(
            self.ProteinSequence, AAP=AAP, AAPName=AAPName
        )
        return res

    def GetMoranAutop(
        self, AAP: Optional[Dict[Any, Any]] = None, AAPName: str = "p"
    ) -> Dict[Any, Any]:
        """
        Moran autocorrelation descriptors for the given property (30).

        Parameters
        ----------
        AAP : Dict[Any, Any]
            contains physicochemical properities of 20 amino acids

        Examples
        --------
        >>> from propy.GetProteinFromUniprot import GetProteinSequence
        >>> protein = GetProteinSequence(ProteinID="Q9NQ39")
        >>> result = GetProDes(protein).GetMoranAutop(AAP={}, AAPName='p')
        """
        if AAP is None:
            AAP = {}
        res = CalculateEachMoranAuto(self.ProteinSequence, AAP=AAP, AAPName=AAPName)
        return res

    def GetGearyAutop(
        self, AAP: Optional[Dict[Any, Any]] = None, AAPName: str = "p"
    ) -> Dict[Any, Any]:
        """
        Geary autocorrelation descriptors for the given property (30).

        Parameters
        ----------
        AAP : Dict[Any, Any]
            contains physicochemical properities of 20 amino acids

        Examples
        --------
        >>> from propy.GetProteinFromUniprot import GetProteinSequence
        >>> protein = GetProteinSequence(ProteinID="Q9NQ39")
        >>> result = GetProDes(protein).GetGearyAutop(AAP={}, AAPName='p')
        """
        if AAP is None:
            AAP = {}
        res = CalculateEachGearyAuto(self.ProteinSequence, AAP=AAP, AAPName=AAPName)
        return res

    def GetSubSeq(self, ToAA: str = "S", window: int = 3) -> List[str]:
        """
        Obtain the sub sequences wit length 2*window+1, whose central point is ToAA.

        ToAA is the central (query point) amino acid in the sub-sequence.

        window is the span.

        Examples
        --------
        >>> from propy.GetProteinFromUniprot import GetProteinSequence
        >>> protein = GetProteinSequence(ProteinID="Q9NQ39")
        >>> result = GetProDes(protein).GetSubSeq(ToAA='S', window=3)
        """
        res = GetSubSequence(self.ProteinSequence, ToAA=ToAA, window=window)
        return res

    def GetALL(
        self,
        paac_lamda: int = 10,
        paac_weight: float = 0.05,
        apaac_lamda: int = 10,
        apaac_weight: float = 0.5,
        socn_maxlag: int = 45,
        qso_maxlag: int = 30,
        qso_weight: float = 0.1,
    ) -> Dict[Any, Any]:
        """
        Calcualte all descriptors except tri-peptide descriptors.

        Parameters
        ----------
        paac_lamda : int, optional (default: 10)
            used by GetPAAC()
            reflects the rank of correlation and is a non-Negative integer,
            such as 15. Note that (1)lamda should NOT be larger than the length
            of input protein sequence; (2) lamda must be non-Negative integer,
            such as 0, 1, 2, ...; (3) when lamda =0, the output of PseAA server
            is the 20-D amino acid composition.
        paac_weight : float, optional (default: 0.05)
            used by GetPAAC()
            is designed for the users to put weight on the additional PseAA
            components with respect to the conventional AA components. The user
            can select any value within the region from 0.05 to 0.7 for the
            weight factor.
        apaac_lamda : int, optional (default: 10)
            Same as "paac_lambda" but for APAAC()
        apaac_weight : float, optional (default: 0.5)
            Same as "paac_weight" but for APAAC()
        socn_maxlag : int, optional  (default: 45)
            Used by GetSOCN()
            is the maximum lag and the length of the protein should be larger
            than maxlag.
        qso_maxlag : int, optional  (default: 30)
            Used by GetQSO()
            is the maximum lag and the length of the protein should be larger
            than maxlag.
        qso_weight : float, optional (default: 0.1)
            Used by GetQSO()
        """

        res: Dict[Any, Any] = {}
        res.update(self.GetAAComp())
        res.update(self.GetDPComp())
        # res.update(self.GetTPComp())
        res.update(self.GetMoreauBrotoAuto())
        res.update(self.GetMoranAuto())
        res.update(self.GetGearyAuto())
        res.update(self.GetCTD())
        res.update(self.GetPAAC(lamda=paac_lamda, weight=paac_weight))
        res.update(self.GetAPAAC(lamda=apaac_lamda, weight=apaac_weight))
        res.update(self.GetSOCN(maxlag=socn_maxlag))
        res.update(self.GetQSO(maxlag=qso_maxlag, weight=qso_weight))
        return res

    def GetAAindex1(self, name: str, path: Optional[str] = ".") -> Dict[str, float]:
        """
        Get the amino acid property values from aaindex1.

        Parameters
        ----------
        name : str
            is the name of amino acid property (e.g., KRIW790103)

        Returns
        -------
        result is a dict form containing the properties of 20 amino acids

        Examples
        --------
        >>> from propy.GetProteinFromUniprot import GetProteinSequence
        >>> protein = GetProteinSequence(ProteinID="Q9NQ39")
        >>> result = GetProDes(protein).GetAAindex1(name="KRIW790103")
        """
        return GetAAIndex1(name, path=path)

    def GetAAindex23(self, name: str, path: Optional[str] = ".") -> Dict[str, float]:
        """
        Get the amino acid property values from aaindex2 and aaindex3.

        Parameters
        ----------
        name is the name of amino acid property (e.g. TANS760101, GRAR740104)

        Returns
        -------
        result is a dict form containing the properties of 400 amino acid pairs

        Examples
        --------
        >>> from propy.GetProteinFromUniprot import GetProteinSequence
        >>> protein = GetProteinSequence(ProteinID="Q9NQ39")
        >>> result = GetProDes(protein).GetAAindex23(name="KRIW790103")
        """
        return GetAAIndex23(name, path=path)
