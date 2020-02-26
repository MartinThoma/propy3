# -*- coding: utf-8 -*-
"""
######################################################################
This is used for testing all commonly used functions in modules.

Authors: Dongsheng Cao and Yizeng Liang.

Date: 2012.09.04

Email: oriental-cds@163.com
######################################################################
"""

modulelists=['AAComposition','Autocorrelation','CTD','QuasiSequenceOrder','PseudoAAC','GetProteinFromUniprot','GetSubSeq']

modules=list()
for i in modulelists:
	modules.append(__import__(i))

AAC  = modules[0]
AC   = modules[1]
CTD  = modules[2]
QSO  = modules[3]
PAAC = modules[4]
GPFU = modules[5]
GSS  = modules[6]

print '...............................................................'

print "testing the GetProteinFromUniprot module"

ProteinSequence = GPFU.GetProteinSequence('P08172')

print '...............................................................'

print "testing the GetSubSeq module"

temp=GSS.GetSubSequence(ProteinSequence,ToAA='D', window=5)

print temp

print '...............................................................'

print "testing the AAComposition module"

temp=AAC.CalculateAAComposition(ProteinSequence)

print temp

temp=AAC.CalculateDipeptideComposition(ProteinSequence)

temp=AAC.GetSpectrumDict(ProteinSequence)

temp=AAC.CalculateAADipeptideComposition(ProteinSequence)

print '...............................................................'

print "testing the Autocorrelation module"


temp=AC.CalculateNormalizedMoreauBrotoAuto(ProteinSequence,[AC._ResidueASA],['ResidueASA'])

print temp

temp=AC.CalculateMoranAuto(ProteinSequence,[AC._ResidueASA],['ResidueASA'])

print temp

temp=AC.CalculateGearyAuto(ProteinSequence,[AC._ResidueASA],['ResidueASA'])

print temp

temp=AC.CalculateAutoTotal(ProteinSequence)

print '...............................................................'

print "testing the CTD module"

temp = CTD.CalculateC(ProteinSequence)

print temp

temp = CTD.CalculateT(ProteinSequence)

print temp

temp = CTD.CalculateD(ProteinSequence)

print temp

temp = CTD.CalculateCTD(ProteinSequence)

print temp

print '...............................................................'

print "testing the QuasiSequenceOrder module"

temp=QSO.GetSequenceOrderCouplingNumberTotal(ProteinSequence,maxlag=30)

print temp

temp=QSO.GetQuasiSequenceOrder(ProteinSequence,maxlag=30,weight=0.1)

print temp

print '...............................................................'

print "testing the PseudoAAC module"

temp= PAAC.GetAPseudoAAC(ProteinSequence,lamda=10,weight=0.5)

print temp


temp= PAAC._GetPseudoAAC(ProteinSequence,lamda=10,weight=0.05)

print temp

print '...............................................................'

print "Tested successfully!"


