# propy3, formerly protpy, is a Python package to compute protein descriptors
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
# Core Library
import sys
import warnings
from typing import List

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
