User Guide
==========


Installation
------------

.. code-block:: bash

   pip install propy3


Download proteins from Uniprot
------------------------------

You can get a protein sequence from the Uniprot website by providing a Uniprot ID:


.. code-block:: python

   from propy.GetProteinFromUniprot import GetProteinSequence as gps
   uniprotid = "P48039"
   proseq = gps(uniprotid)

   print(proseq)

gives

.. code-block

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

.. code-block

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

.. code-block

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
------------------------------------------------

You could get the properties of amino acids from the AAindex database by
providing a property name (e.g., KRIW790103). The output is given in the form
of dictionary.

If the user provides the directory containing the AAindex database (the AAindex
database could be downloaded from
ftp://ftp.genome.jp/pub/db/community/aaindex/. It consists of three files:
aaindex1, aaindex2 and aaindex3), the program will read the given database to
get the property.

TODO

Calculating protein descriptors
-------------------------------

TODO
