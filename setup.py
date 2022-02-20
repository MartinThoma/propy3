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
"""Set up the propy package"""

# Third party
from setuptools import setup

packagedata = {
    "propy": [
        "aaindex1",
        "aaindex2",
        "aaindex3",
        "html/*",
        "instruction/*",
        "data/*",
        "aaindex/*",
    ]
}


setup(package_data=packagedata)
