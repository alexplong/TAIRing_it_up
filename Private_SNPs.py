#!/usr/bin/python

import csv
import sys

SNPs = 'ALL.tab'

Bonly = 0
Conly = 0
Lonly = 0
Tonly = 0

BC = 0
BL = 0
BT = 0
CL = 0
CT = 0
LT = 0

BCL = 0
BCT = 0
BLT = 0
CLT = 0

Bref = 0
Cref = 0
Lref = 0
Tref = 0

ALL = 0
SHR = 0
TOT = 0

with open(SNPs, 'rb') as SNPfile:
    line = csv.reader(SNPfile, delimiter = '\t')
    for row in line:
        GEN = row[0]
        POS = row[1]
        REF = row[2]
        BUR = row[3]
        CAN = row[4]
        LER = row[5]
        TSU = row[6]

        TOT += 1

        if (BUR[0] == CAN[0] and BUR[0] == LER[0] and BUR[0] == TSU[0] and BUR[0] != REF):
            ALL += 1
        elif (BUR[0] != CAN[0] and BUR[0] != LER[0] and BUR[0] != TSU[0]) and BUR[0] != REF:
            Bonly += 1
        elif (CAN[0] != BUR[0] and CAN[0] != LER[0] and CAN[0] != TSU[0]) and CAN[0] != REF:
            Conly += 1
        elif (LER[0] != CAN[0] and LER[0] != BUR[0] and LER[0] != TSU[0]) and LER[0] != REF:
            Lonly += 1
        elif (TSU[0] != CAN[0] and TSU[0] != LER[0] and TSU[0] != BUR[0]) and TSU[0] != REF:
            Tonly += 1
        elif (BUR[0] == CAN[0] and BUR[0] != LER[0] and BUR[0] != TSU[0]) and BUR[0] != REF:
            BC += 1
        elif (BUR[0] == LER[0] and BUR[0] != CAN[0] and BUR[0] != TSU[0]) and BUR[0] != REF:
            BL += 1
        elif (BUR[0] == TSU[0] and BUR[0] != LER[0] and BUR[0] != CAN[0]) and BUR[0] != REF:
            BT += 1
        elif (LER[0] == CAN[0] and LER[0] != BUR[0] and LER[0] != TSU[0]) and LER[0] != REF:
            CL += 1
        elif (TSU[0] == CAN[0] and TSU[0] != LER[0] and TSU[0] != BUR[0]) and TSU[0] != REF:
            CT += 1
        elif (LER[0] == TSU[0] and LER[0] != BUR[0] and LER[0] != CAN[0]) and LER[0] != REF:
            LT +=1
        elif (BUR[0] == CAN[0] and BUR[0] == LER[0] and BUR[0] != TSU[0]) and BUR[0] != REF:
            BCL += 1
        elif (BUR[0] == CAN[0] and CAN[0] == TSU[0] and BUR[0] != LER[0]) and BUR[0] != REF:
            BCT += 1
        elif (BUR[0] == LER[0] and BUR[0] == TSU[0] and BUR[0] != CAN[0]) and BUR[0] != REF:
            BLT += 1
        elif (CAN[0] == LER[0] and CAN[0] == TSU[0] and CAN[0] != BUR[0]) and CAN[0] != REF:
            CLT += 1
        elif (REF == CAN[0] and REF == LER[0] and REF == TSU[0] and REF == BUR[0]):
            SHR += 1
        if (BUR[0] == REF and CAN[0] != REF and LER[0] != REF and TSU[0] != REF):
            Bref += 1
        if (CAN[0] == REF and BUR[0] != REF and LER[0] != REF and TSU[0] != REF):
            Cref += 1
        if (LER[0] == REF and BUR[0] != REF and CAN[0] != REF and TSU[0] != REF):
            Lref += 1
        if (TSU[0] == REF and BUR[0] != REF and CAN[0] != REF and LER[0] != REF):
            Tref += 1



print 'Total Number of SNPs:', TOT
print '##############################'
print '# of Private SNPs'
print 'Bur_0' '\t' 'Can_0' '\t' 'Ler_0' '\t' 'Tsu_0'
print Bonly, '\t', Conly, '\t', Lonly, '\t', Tonly
print '##############################'
print '# of Shared SNPs (Pairs)'
#print 'Bur_0/Can_0', '\t', 'Bur_0/Ler_0', '\t', 'Bur_0/Tsu_0', '\t', 'Can_0/Ler_0', '\t', 'Can_0/Tsu_0', '\t', 'Ler_0/Tsu_0'
print 'BC', '\t', 'BL', '\t', 'BT', '\t', 'CL', '\t', 'CT', '\t', 'LT'
print BC, '\t', BL, '\t', BT, '\t', CL, '\t', CT, '\t', LT
print '##############################'
print '# of Shared SNPs (Triplets)'
#print 'Bur_0/Can_0/Ler_0', '\t', 'Bur_0/Can_0/Tsu_0', '\t', 'Bur_0/Ler_0/Tsu_0', '\t', 'Can_0/Ler_0/Tsu_0'
print 'BCL', '\t', 'BCT', '\t', 'BLT', '\t', 'CLT'
print BCL, '\t', BCT, '\t', BLT, '\t', CLT
print '##############################'
print '# of SNPs shared with REF'
print 'Bur_0' '\t' 'Can_0' '\t' 'Ler_0' '\t' 'Tsu_0'
print Bref, '\t', Cref, '\t', Lref, '\t', Tref
