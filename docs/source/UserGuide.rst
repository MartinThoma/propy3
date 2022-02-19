Installation
============

.. code-block:: bash

   pip install propy3


Download proteins from Uniprot
==============================

You can get a protein sequence from the Uniprot website by providing a Uniprot ID:


.. code-block:: python

   from propy.GetProteinFromUniprot import GetProteinSequence as gps

   uniprotid = "P48039"
   proseq = gps(uniprotid)

   print(proseq)

gives

.. code-block:: text

   MQGNGSALPNASQPVLRGDGARPSWLASALACVLIFTIVVDILGNLLVILSVYRNKKLRNAGNIFVVSLAVA\
   DLVVAIYPYPLVLMSIFNNGWNLGYLHCQVSGFLMGLSVIGSIFNITGIAINRYCYICHSLKYDKLYSSKNS\
   LCYVLLIWLLTLAAVLPNLRAGTLQYDPRIYSCTFAQSVSSAYTIAVVVFHFLVPMIIVIFCYLRIWILVLQ\
   VRQRVKPDRKPKLKPQDFRNFVTMFVVFVLFAICWAPLNFIGLAVASDPASMVPRIPEWLFVASYYMAYFNS\
   CLNAIIYGLLNQNFRKEYRRIIVSLCTARVFFVDSSNDVADRVKWKPSPLMTNNNVVKVDSV


You can get the window✕ 2 + 1 sub-sequences whose central point is the given
amino acid ToAA.

.. code-block:: python

   from propy import GetSubSeq

   subseq = GetSubSeq.GetSubSequence(proseq, ToAA="S", window=5)
   print(subseq)

gives

.. code-block:: python

   ['MQGNGSALPNA', 'ALPNASQPVLR', 'DGARPSWLASA', 'PSWLASALACV', 'LLVILSVYRNK',
    'NIFVVSLAVAD', 'PLVLMSIFNNG', 'LHCQVSGFLMG', 'FLMGLSVIGSI', 'LSVIGSIFNIT',
    'CYICHSLKYDK', 'YDKLYSSKNSL', 'DKLYSSKNSLC', 'YSSKNSLCYVL', 'DPRIYSCTFAQ',
    'CTFAQSVSSAY', 'FAQSVSSAYTI', 'AQSVSSAYTIA', 'GLAVASDPASM', 'ASDPASMVPRI',
    'WLFVASYYMAY', 'MAYFNSCLNAI', 'RRIIVSLCTAR', 'VFFVDSSNDVA', 'FFVDSSNDVAD',
    'VKWKPSPLMTN']

You can also get several protein sequences by providing a file containing
Uniprot IDs of these proteins.

.. code-block:: python

   from propy.GetProteinFromUniprot import GetProteinSequenceFromTxt as gpst

   tag = gpst("propy/data", "target.txt", "target1.txt")

prints

.. code-block:: text

    --------------------------------------------------------------------------------
    The 1 protein sequence has been downloaded!
    MADSCRNLTYVRGSVGPATSTLMFVAGVVGNGLALGILSARRPARPSAFAVLVTGLAATDLLGTSFLSPAVFVAYARNSSLLGLARGGPALCDAFAFAMTFFGLASMLILFAMAVERCLALSHPYLYAQLDGPRCARLALPAIYAFCVLFCALPLLGLGQHQQYCPGSWCFLRMRWAQPGGAAFSLAYAGLVALLVAAIFLCNGSVTLSLCRMYRQQKRHQGSLGPRPRTGEDEVDHLILLALMTVVMAVCSLPLTIRCFTQAVAPDSSSEMGDLLAFRFYAFNPILDPWVFILFRKAVFQRLKLWVCCLCLGPAHGDSQTPLSQLASGRRDPRAPSAPVGKEGSCVPLSAWGEGQVEPLPPTQQSSGSAVGTSSKAEASVACSLC
    --------------------------------------------------------------------------------

    TODO: HTTP Error 300!

The downloaded protein sequences have been saved in "propy/data/target1.txt".

You could check whether the input sequence is a valid protein sequence or not:

.. code-block:: python

   from propy import ProCheck

   temp = ProCheck.ProteinCheck(proseq)
   print(tmp)

which prints :code:`350`. This output is the number of the protein sequence if
it is valid; otherwise 0.


Obtaining the property from the AAindex database
================================================

You could get the properties of amino acids from the AAindex database by
providing a property name (e.g., KRIW790103). The output is given in the form
of dictionary.

If the user provides the directory containing the AAindex database (the AAindex
database could be downloaded from
ftp://ftp.genome.jp/pub/db/community/aaindex/. It consists of three files:
aaindex1, aaindex2 and aaindex3), the program will read the given database to
get the property.

.. code-block:: pycon
   >>> from propy.AAIndex import GetAAIndex1, GetAAIndex23
   >>> temp1 = GetAAIndex1("KRIW790103")
   >>> temp1
   {'A': 27.5, 'R': 105.0, 'N': 58.7, 'D': 40.0, 'C': 44.6, 'E': 62.0, 'Q': 80.7, 'G': 0.0, 'H': 79.0, 'I': 93.5, 'L': 93.5, 'K': 100.0, 'M': 94.1, 'F': 115.5, 'P': 41.9, 'S': 29.3, 'T': 51.3, 'W': 145.5, 'Y': 117.3, 'V': 71.5}



Calculating protein descriptors
===============================

.. code-block:: pycon
   >>> from propy import PyPro
   >>> from propy.GetProteinFromUniprot import GetProteinSequence as gps
   >>> proteinsequence = gps("P48039")
   >>> DesObject = PyPro.GetProDes(proteinsequence)
   >>> print(DesObject.GetCTD())
   {'_PolarizabilityC1': 0.257, '_PolarizabilityC2': 0.494, '_PolarizabilityC3': 0.249, '_SolventAccessibilityC1': 0.546, '_SolventAccessibilityC2': 0.22, '_SolventAccessibilityC3': 0.234, '_SecondaryStrC1': 0.351, '_SecondaryStrC2': 0.38, '_SecondaryStrC3': 0.269, '_ChargeC1': 0.089, '_ChargeC2': 0.871, '_ChargeC3': 0.04, '_PolarityC1': 0.5, '_PolarityC2': 0.271, '_PolarityC3': 0.229, '_NormalizedVDWVC1': 0.334, '_NormalizedVDWVC2': 0.417, '_NormalizedVDWVC3': 0.249, '_HydrophobicityC1': 0.22, '_HydrophobicityC2': 0.331, '_HydrophobicityC3': 0.449, '_PolarizabilityT12': 0.258, '_PolarizabilityT13': 0.097, '_PolarizabilityT23': 0.289, '_SolventAccessibilityT12': 0.209, '_SolventAccessibilityT13': 0.246, '_SolventAccessibilityT23': 0.1, '_SecondaryStrT12': 0.261, '_SecondaryStrT13': 0.212, '_SecondaryStrT23': 0.149, '_ChargeT12': 0.143, '_ChargeT13': 0.011, '_ChargeT23': 0.069, '_PolarityT12': 0.266, '_PolarityT13': 0.198, '_PolarityT23': 0.129, '_NormalizedVDWVT12': 0.281, '_NormalizedVDWVT13': 0.155, '_NormalizedVDWVT23': 0.232, '_HydrophobicityT12': 0.149, '_HydrophobicityT13': 0.16, '_HydrophobicityT23': 0.292, '_PolarizabilityD1001': 0.857, '_PolarizabilityD1025': 19.429, '_PolarizabilityD1050': 44.857, '_PolarizabilityD1075': 74.857, '_PolarizabilityD1100': 99.714, '_PolarizabilityD2001': 0.571, '_PolarizabilityD2025': 23.429, '_PolarizabilityD2050': 48.0, '_PolarizabilityD2075': 72.571, '_PolarizabilityD2100': 100.0, '_PolarizabilityD3001': 0.286, '_PolarizabilityD3025': 33.143, '_PolarizabilityD3050': 58.571, '_PolarizabilityD3075': 78.571, '_PolarizabilityD3100': 98.857, '_SolventAccessibilityD1001': 0.857, '_SolventAccessibilityD1025': 21.714, '_SolventAccessibilityD1050': 45.714, '_SolventAccessibilityD1075': 71.429, '_SolventAccessibilityD1100': 100.0, '_SolventAccessibilityD2001': 0.571, '_SolventAccessibilityD2025': 26.0, '_SolventAccessibilityD2050': 62.286, '_SolventAccessibilityD2075': 85.429, '_SolventAccessibilityD2100': 99.429, '_SolventAccessibilityD3001': 0.286, '_SolventAccessibilityD3025': 28.286, '_SolventAccessibilityD3050': 50.0, '_SolventAccessibilityD3075': 76.286, '_SolventAccessibilityD3100': 99.714, '_SecondaryStrD1001': 0.286, '_SecondaryStrD1025': 21.143, '_SecondaryStrD1050': 47.143, '_SecondaryStrD1075': 72.571, '_SecondaryStrD1100': 98.857, '_SecondaryStrD2001': 4.286, '_SecondaryStrD2025': 30.0, '_SecondaryStrD2050': 54.0, '_SecondaryStrD2075': 73.143, '_SecondaryStrD2100': 100.0, '_SecondaryStrD3001': 0.857, '_SecondaryStrD3025': 19.429, '_SecondaryStrD3050': 41.143, '_SecondaryStrD3075': 76.286, '_SecondaryStrD3100': 99.714, '_ChargeD1001': 4.857, '_ChargeD1025': 35.714, '_ChargeD1050': 62.857, '_ChargeD1075': 86.571, '_ChargeD1100': 98.857, '_ChargeD2001': 0.286, '_ChargeD2025': 24.286, '_ChargeD2050': 47.714, '_ChargeD2075': 73.143, '_ChargeD2100': 100.0, '_ChargeD3001': 5.429, '_ChargeD3025': 20.857, '_ChargeD3050': 66.571, '_ChargeD3075': 87.143, '_ChargeD3100': 99.429, '_PolarityD1001': 0.286, '_PolarityD1025': 25.429, '_PolarityD1050': 50.0, '_PolarityD1075': 71.714, '_PolarityD1100': 100.0, '_PolarityD2001': 0.857, '_PolarityD2025': 19.429, '_PolarityD2050': 45.143, '_PolarityD2075': 74.286, '_PolarityD2100': 99.714, '_PolarityD3001': 0.571, '_PolarityD3025': 26.857, '_PolarityD3050': 61.714, '_PolarityD3075': 85.429, '_PolarityD3100': 99.429, '_NormalizedVDWVD1001': 0.857, '_NormalizedVDWVD1025': 20.857, '_NormalizedVDWVD1050': 47.143, '_NormalizedVDWVD1075': 74.857, '_NormalizedVDWVD1100': 99.714, '_NormalizedVDWVD2001': 0.571, '_NormalizedVDWVD2025': 21.714, '_NormalizedVDWVD2050': 46.286, '_NormalizedVDWVD2075': 72.571, '_NormalizedVDWVD2100': 100.0, '_NormalizedVDWVD3001': 0.286, '_NormalizedVDWVD3025': 33.143, '_NormalizedVDWVD3050': 58.571, '_NormalizedVDWVD3075': 78.571, '_NormalizedVDWVD3100': 98.857, '_HydrophobicityD1001': 0.571, '_HydrophobicityD1025': 26.0, '_HydrophobicityD1050': 62.286, '_HydrophobicityD1075': 85.429, '_HydrophobicityD1100': 99.429, '_HydrophobicityD2001': 0.857, '_HydrophobicityD2025': 22.857, '_HydrophobicityD2050': 45.143, '_HydrophobicityD2075': 74.286, '_HydrophobicityD2100': 99.714, '_HydrophobicityD3001': 0.286, '_HydrophobicityD3025': 25.143, '_HydrophobicityD3050': 51.143, '_HydrophobicityD3075': 71.429, '_HydrophobicityD3100': 100.0}
   >>> print(DesObject.GetAAComp())
   {'A': 7.714, 'R': 5.143, 'N': 6.571, 'D': 3.429, 'C': 2.857, 'E': 0.571, 'Q': 2.571, 'G': 4.286, 'H': 0.857, 'I': 8.0, 'L': 12.286, 'K': 3.714, 'M': 2.286, 'F': 5.714, 'P': 4.857, 'S': 7.714, 'T': 2.571, 'W': 2.0, 'Y': 5.143, 'V': 11.714}
   >>> paac = DesObject.GetPAAC(
   ...         lamda=10, weight=0.05
   ...     )
   >>> 
   >>> for i in paac:
   ...         print(i)
   ... 
   PAAC1
   PAAC2
   PAAC3
   PAAC4
   PAAC5
   PAAC6
   PAAC7
   PAAC8
   PAAC9
   PAAC10
   PAAC11
   PAAC12
   PAAC13
   PAAC14
   PAAC15
   PAAC16
   PAAC17
   PAAC18
   PAAC19
   PAAC20
   PAAC21
   PAAC22
   PAAC23
   PAAC24
   PAAC25
   PAAC26
   PAAC27
   PAAC28
   PAAC29
   PAAC30
