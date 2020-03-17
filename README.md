[![PyPI version](https://badge.fury.io/py/propy3.svg)](https://badge.fury.io/py/propy3)
[![Python Support](https://img.shields.io/pypi/pyversions/propy3.svg)](https://pypi.org/project/propy3/)
[![Documentation Status](https://readthedocs.org/projects/propy3/badge/?version=latest)](https://propy3.readthedocs.io/en/latest/?badge=latest)
[![Build Status](https://travis-ci.org/MartinThoma/propy3.svg?branch=master)](https://travis-ci.org/MartinThoma/propy3)
[![Coverage Status](https://coveralls.io/repos/github/MartinThoma/propy3/badge.svg?branch=master)](https://coveralls.io/github/MartinThoma/propy3?branch=master)

# propy3

`propy3` is a drop-in replacement for [`propy`](https://code.google.com/archive/p/protpy/).
The original project was developed by Dongsheng Cao. See the commit history
for all changes made afterwards.

The reason for creating this fork of propy is to add Python 3 support.

The only point where you have to enter `propy3` is at installation. Afterwards,
you simply `import propy`.

## Introduction

Sequence-derived structural and physicochemical features are highly useful for
representing and distinguishing proteins or peptides of different structural,
functional and interaction properties, and have been extensively used in
developing methods and software for predicting protein structural and
functional classes, protein-protein interactions, drug-target interactions,
protein substrates, molecualr binding sites on proteins, subcellular locations,
protein crystallization propensity and peptides of specific properties. In
order to conveniently apply these structural features from a protein sequence
for researchers, we developed a propy package using pure python language, which
could calculate a large number of protein descriptors from a protein sequence.

## Features

The propy package has the following significant features:

1. It is written by the pure python language. It only needs the support of some
   built-in modules in the python software.
2. For academic users, it is free of charge. They can freely use and distribute
   it. For commercial purpose, they must contact the author.
3. It can calculate a large number of protein descriptors including: amino acid
   composition descriptors, dipeptide composition descriptors, tri-peptide
   composition descriptors, Normalized Moreau-Broto autocorrelation
   descriptors, Moran autocorrelation descriptors, Geary autocorrelation
   descriptors, Composition, Transition, Distribution descriptors (CTD),
   sequence order coupling numbers, quasi-sequence order descriptors, pseudo
   amino acid composition descriptors, amphiphilic pseudo amino acid
   composition descriptors.
4. The users could specify the needed properties of 20 amino acids to calculate
   the corresponding protein descriptors.
5. The package includes the module which could directly download the protein
   sequence form uniprot website by uniprot id.
6. The package includes the module which could automatrically download the
   property from the AAindex database. Thus, the user could calcualte thousands
   of protein features.

The protein descriptors calculated by propy

1. AAC: amino acid composition descriptors (20)
2. DPC: dipeptide composition descriptors (400)
3. TPC: tri-peptide composition descriptors (8000)
4. MBauto: Normalized Moreau-Broto autocorrelation descriptors (depend on the given properties, the default is 240)
5. Moranauto: Moran autocorrelation descriptors(depend on the given properties, the default is 240)
6. Gearyauto: Geary autocorrelation descriptors(depend on the given properties, the default is 240)
6. CTD: Composition, Transition, Distribution descriptors (CTD) (21+21+105=147)
7. SOCN: sequence order coupling numbers (depend on the choice of maxlag, the default is 60)
8. QSO: quasi-sequence order descriptors (depend on the choice of maxlag, the default is 100)
9. PAAC: pseudo amino acid composition descriptors (depend on the choice of lamda, the default is 50)
10. APAAC: amphiphilic pseudo amino acid composition descriptors(depend on the choice of lamda, the default is 50)

## Install

```
pip install propy3
```

## Usage Example

For more examples, please see the user guide.

```python
from propy import PyPro
from propy.GetProteinFromUniprot import GetProteinSequence

# download the protein sequence by uniprot id
proteinsequence = GetProteinSequence("P48039")

DesObject = PyPro.GetProDes(proteinsequence)  # construct a GetProDes object
print(DesObject.GetCTD())  # calculate 147 CTD descriptors
print(DesObject.GetAAComp())  # calculate 20 amino acid composition descriptors

# calculate 30 pseudo amino acid composition descriptors
paac = DesObject.GetPAAC(lamda=10, weight=0.05)

for i in paac:
    print(i, paaci)
```
