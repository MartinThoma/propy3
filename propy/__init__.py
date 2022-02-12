# Core Library
from typing import List
import sys
import warnings

_python_version = sys.version_info

if _python_version.major == 2:
    warnings.warn("Python 2 is not supported. Please use Python 3.")
if _python_version.major == 3 and _python_version.minor < 8:
    warnings.warn(
        "Python 3.6 and Python 3.7 might get deprecated. "
        "Please participate in the discussion: "
        "https://github.com/MartinThoma/propy3/issues/12"
    )

AALetter: List[str] = list("ARNDCEQGHILKMFPSTWYV")

ProteinSequence_docstring = """ProteinSequence: str
        a pure protein sequence"""
