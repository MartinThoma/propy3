# -*- coding: utf-8 -*-
INSTRUCTION
Sequence-derived structural and physicochemical features are highly useful for representing 
and distinguishing proteins or peptides of different structural, functional and interaction 
properties, and have been extensively used in developing methods and software for predicting
protein structural and functional classes, protein-protein interactions, drug-target interactions,
protein substrates, molecualr binding sites on proteins, subcellular locations, protein crystallization 
propensity and peptides of specific properties. In order to conveniently apply these structural features
from a protein sequence for researchers, we developed a propy package using pure python language, which 
could calculate a large number of protein descriptors from a protein sequence.

#############################################FEATURES##################################################

The propy package has the following significant features:
(1): It is written by the pure python language. It only needs the support 
of some built-in modules in the python software.
(2): For academic users, it is free of charge. They can freely use and distribute
 it. For commercial purpose, they must contact the author.
(3): It can calculate a large number of protein descriptors including: 
amino acid composition descriptors, dipeptide composition descriptors, 
tri-peptide composition descriptors, Normalized Moreau-Broto autocorrelation descriptors, 
Moran autocorrelation descriptors, Geary autocorrelation descriptors, Composition, 
Transition, Distribution descriptors (CTD), sequence order coupling numbers, 
quasi-sequence order descriptors, pseudo amino acid composition descriptors, 
amphiphilic pseudo amino acid composition descriptors.
(4): The users could specify the needed properties of 20 amino acids to calculate 
the corresponding protein descriptors.
(5): The package includes the module which could directly download the protein 
sequence form uniprot website by uniprot id. 
(6): The package includes the module which could automatrically download the property
from the AAindex database. Thus, the user could calcualte thousands of protein features.

##########################################################################################################
The protein descriptors calculated by propy 

(1) AAC: amino acid composition descriptors (20)
(2) DPC: dipeptide composition descriptors (400)
(3) TPC: tri-peptide composition descriptors (8000)
(4) MBauto: Normalized Moreau-Broto autocorrelation descriptors (depend on the given properties, the default is 240)
(5) Moranauto: Moran autocorrelation descriptors(depend on the given properties, the default is 240)
(6) Gearyauto: Geary autocorrelation descriptors(depend on the given properties, the default is 240)
(6) CTD: Composition, Transition, Distribution descriptors (CTD) (21+21+105=147)
(7) SOCN: sequence order coupling numbers (depend on the choice of maxlag, the default is 60)
(8) QSO: quasi-sequence order descriptors (depend on the choice of maxlag, the default is 100)
(9) PAAC: pseudo amino acid composition descriptors (depend on the choice of lamda, the default is 50)
(10) APAAC: amphiphilic pseudo amino acid composition descriptors(depend on the choice of lamda, the default is 50) 
##########################################################################################################
Download

propy can be download from http://protpy.googlecode.com/files/propy-1.0.tar.gz
##########################################################################################################
Install

On Windows:

(1): download the propy package (.gz)

(2): extract or uncompress the .gz file

(3): cd propy-1.0

(4): python setup.py install

On Linux:

(1): download the propy package (.tar.gz)

(2): tar -zxf propy-1.0.tar.gz

(3): cd propy-1.0

(4): python setup.py install or sudo python setup.py install

##########################################################################################################

Example:

For more examples, please see the user guide.

from propy import PyPro

from propy.GetProteinFromUniprot import GetProteinSequence

proteinsequence=GetProteinSequence('P33765')      ##download the protein sequence by uniprot id

DesObject=PyPro.GetProDes(proteinsequence)        ##construct a GetProDes object

print DesObject.GetCTD()                          ##calculate 147 CTD descriptors

print DesObject.GetAAComp()                       ##calculate 20 amino acid composition descriptors

paac=DesObject.GetPAAC(lamda=10,weight=0.05)      ##calculate 30 pseudo amino acid composition descriptors

for i in paac:

    print i, paaci 

##########################################################################################################
