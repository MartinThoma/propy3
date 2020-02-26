# -*- coding: utf-8 -*-
"""
Set up the propy package

Authors: Dongsheng Cao and Yizeng Liang.

Date: 2012.09.11

Email: oriental-cds@163.com
"""

from distutils.core import setup 


#datafiles=[('propy/html',['src/html/AAComposition.html','src/html/Autocorrelation.html','src/html/CTD.html','src/html/GetSubSeq.html','src/html/GetProteinFromUniprot.html','src/html/PseudoAAC.html','src/html/QuasiSequenceOrder.html']),('propy',['README.txt']),('propy/docs',['src/instruction/Schneider-Wrede distance.xls','src/instruction/UserGuide.pdf','src/instruction/Manual.pdf','src/instruction/Grantham.xls']),('propy/data',['src/data/target.txt']),('propy/aaindexa',['src/aaindex/aaindex1','src/aaindex/aaindex2','src/aaindex/aaindex3','src/aaindex/changelog','src/aaindex/aaindex.doc','src/aaindex/Fig.4.GIF','src/aaindex/Fig.5-1.GIF','src/aaindex/Fig.5-2.GIF']),('',['src/propy/aaindex1','src/propy/aaindex2','src/propy/aaindex3'])]


packagedata={'propy': ['aaindex1','aaindex2','aaindex3','html/*','instruction/*','data/*','aaindex/*']}


setup(name = 'propy', 

	version = '1.0', 
	
	description ="Compute protein descriptors",
	
	author = "Dongsheng Cao",
	
	author_email = "oriental-cds@163.com",
	
	url ="http://cbdd.csu.edu.cn/index",
	
	license = "GPL",
	
	packages = ['propy'],
	
	package_data=packagedata,
	
#	data_files = datafiles,
	
	package_dir={'propy':'src/propy'},
	
	scripts = [],
	
	py_modules = []

	)

