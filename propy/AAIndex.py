"""
This module is used for obtaining the properties of amino acids or their pairs
from the aaindex database.

Authors: Dongsheng Cao and Yizeng Liang.
Date: 2012.09.10
Email: oriental-cds@163.com
"""

# Core Library
import logging
import os
import sys
from typing import Any, Dict, List, Optional

# Third party
import pkg_resources

logger = logging.getLogger(__name__)

AALetter: List[str] = list("ARNDCEQGHILKMFPSTWYV")

_aaindex: Dict[Any, Any] = {}


class Record:
    """Amino acid index (AAindex) Record."""

    aakeys = "ARNDCQEGHILKMFPSTWYV"

    def __init__(self):
        self.key = None
        self.desc = ""
        self.ref = ""
        self.authors = ""
        self.title = ""
        self.journal = ""
        self.correlated = {}
        self.index = {}
        self.comment = ""

    def extend(self, row):
        i = len(self.index)
        for x in row:
            self.index[self.aakeys[i]] = x
            i += 1

    def get(self, aai, aaj=None, d=None):
        assert aaj is None
        return self.index.get(aai, d)

    def __getitem__(self, aai):
        return self.get(aai)

    def median(self):
        x = sorted([_f for _f in list(self.index.values()) if _f])
        half = len(x) // 2
        if len(x) % 2 == 1:
            return x[half]
        return (x[half - 1] + x[half]) / 2.0

    def __str__(self):
        desc = self.desc.replace("\n", " ").strip()
        return "%s(%s: %s)" % (self.__class__.__name__, self.key, desc)


class MatrixRecord(Record):
    """Matrix record for mutation matrices or pair-wise contact potentials."""

    def __init__(self):
        Record.__init__(self)
        self.index: List[Any] = []  # type: ignore
        self.rows = {}
        self.cols = {}

    def extend(self, row):
        self.index.append(row)

    def _get(self, aai, aaj):
        i = self.rows[aai]
        j = self.cols[aaj]
        return self.index[i][j]

    def get(self, aai, aaj, d=None):
        try:
            return self._get(aai, aaj)
        except Exception as e:
            logger.debug(e)
            pass
        try:
            return self._get(aaj, aai)
        except Exception as e:
            logger.debug(e)
            return d

    def __getitem__(self, aaij):
        return self.get(aaij[0], aaij[1])

    def median(self):
        x = []
        for y in self.index:
            x.extend([_f for _f in y if _f])
        x.sort()
        if len(x) % 2 == 1:
            return x[len(x) // 2]
        return sum(x[len(x) // 2 - 1 : len(x) // 2 + 1]) / 2.0


def search(pattern, searchtitle=True, casesensitive=False):
    """
    Search for pattern in description and title (optional) of all records and
    return matched records as list. By default search case insensitive.
    """
    whatcase = lambda i: i
    if not casesensitive:
        pattern = pattern.lower()
        whatcase = lambda i: i.lower()
    matches = []
    for record in _aaindex.values():
        if (
            pattern in whatcase(record.desc)
            or searchtitle
            and pattern in whatcase(record.title)
        ):
            matches.append(record)
    return matches


def grep(pattern):
    """
    Search for pattern in title and description of all records (case
    insensitive) and print results on standard output.
    """
    for record in search(pattern):
        print(record)


def get(key: str):
    """Get record for key."""
    if len(_aaindex) == 0:
        init()
    return _aaindex[key]


def _float_or_None(x) -> Optional[float]:
    if x == "NA" or x == "-":
        return None
    return float(x)


def init(path=None, index="123"):
    """
    Read in the aaindex files. You need to run this (once) before you can
    access any records. If the files are not within the current directory, you
    need to specify the correct directory path. By default all three aaindex
    files are read in.
    """
    index = str(index)
    if path is None:
        filepath = pkg_resources.resource_filename(__name__, "aaindex1")
        path = os.path.dirname(filepath)
        print("path =", path, file=sys.stderr)
    if "1" in index:
        _parse(os.path.join(path, "aaindex1"), Record)
    if "2" in index:
        _parse(os.path.join(path, "aaindex2"), MatrixRecord)
    if "3" in index:
        _parse(os.path.join(path, "aaindex3"), MatrixRecord)


def init_from_file(filename, type=Record):
    _parse(filename, type)


def _parse(filename: str, rec, quiet: bool = True):
    """
    Parse aaindex input file. `rec` must be `Record` for aaindex1 and
    `MarixRecord` for aaindex2 and aaindex3.
    """
    if not os.path.exists(filename):
        import urllib.request
        import urllib.parse
        import urllib.error

        url = (
            "ftp://ftp.genome.jp/pub/db/community/aaindex/" + os.path.split(filename)[1]
        )
        logger.debug(f'Downloading "{url}"')
        filename = urllib.request.urlretrieve(url, filename)[0]
    logger.debug(f'Saved to "{filename}"')
    f = open(filename)

    current = rec()
    lastkey = None
    for line in f:
        key = line[0:2]
        if key[0] == " ":
            key = lastkey  # type: ignore
        if key == "//":
            _aaindex[current.key] = current
            current = rec()
        elif key == "H ":
            current.key = line[2:].strip()
        elif key == "R ":
            current.ref += line[2:]
        elif key == "D ":
            current.desc += line[2:]
        elif key == "A ":
            current.authors += line[2:]
        elif key == "T ":
            current.title += line[2:]
        elif key == "J ":
            current.journal += line[2:]
        elif key == "* ":
            current.comment += line[2:]
        elif key == "C ":
            a = line[2:].split()
            for i in range(0, len(a), 2):
                current.correlated[a[i]] = float(a[i + 1])
        elif key == "I ":
            a = line[1:].split()
            if a[0] != "A/L":
                current.extend([_float_or_None(el) for el in a])
            elif list(Record.aakeys) != [i[0] for i in a] + [i[-1] for i in a]:
                print("Warning: wrong amino acid sequence for", current.key)
            else:
                try:
                    assert list(Record.aakeys[:10]) == [i[0] for i in a]
                    assert list(Record.aakeys[10:]) == [i[2] for i in a]
                except Exception as e:
                    logger.debug(e)
                    print("Warning: wrong amino acid sequence for", current.key)
        elif key == "M ":
            a = line[2:].split()
            if a[0] == "rows":
                if a[4] == "rows":
                    a.pop(4)
                assert a[3] == "cols" and len(a) == 6
                i = 0
                for aa in a[2]:
                    current.rows[aa] = i
                    i += 1
                i = 0
                for aa in a[5]:
                    current.cols[aa] = i
                    i += 1
            else:
                current.extend([_float_or_None(el) for el in a])
        elif not quiet:
            print('Warning: line starts with "%s"' % (key))
        lastkey = key
    f.close()


def GetAAIndex1(name: str, path: str = ".") -> Dict[str, float]:
    """
    Get the amino acid property values from aaindex1.

    Parameters
    ----------
    name : str
        name of amino acid property (e.g., KRIW790103)

    Returns
    -------
    result : Dict[str, float]
        contains the properties of 20 amino acids

    Examples
    --------
    >>> result = GetAAIndex1("KRIW790103")
    """

    init(path=path)
    name = str(name)
    temp = get(name.strip())
    res = {}
    for i in AALetter:
        res[i] = temp.get(i)
    return res


def GetAAIndex23(name: str, path: str = ".") -> Dict[str, float]:
    """
    Get the amino acid property values from aaindex2 and aaindex3.

    Parameters
    ----------
    name : str
        name of amino acid property (e.g.,TANS760101,GRAR740104)

    Returns
    -------
    result : Dict[str, float]
        contains the properties of 400 amino acid pairs

    Examples
    --------
    >>> result = GetAAIndex23("TANS760101")
    """
    init(path=path)
    name = str(name)
    temp = get(name.strip())
    res = {}
    for i in AALetter:
        for j in AALetter:
            res[i + j] = temp.get(i, j)
    return res
