# -*- coding: utf-8 -*-
"""
Set up the propy package

Authors: Dongsheng Cao and Yizeng Liang.

Date: 2012.09.11

Email: oriental-cds@163.com
"""

# Core Library
from distutils.core import setup

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
