from numpy import *
import subprocess
import os, sys
import myfunc.util as util

#*** Read Caselist **
listPath = '/home/utsumi/bin/PMM/ret-epc/caselist.csv'
f=open(listPath); lines=f.readlines(); f.close()
lcase = []
for line in lines[1:]:
    line = line.strip().split(',')
    Year,Mon,Day=map(int, line[:3])
    clat = float(line[3])
    clon = float(line[4])
    oid  = int(line[5])
    lcase.append([Year,Mon,Day,clat,clon,oid])
#********************
#lcase = lcase[:1]  # test
for case in lcase:
    print case
    Year,Mon,Day,clat,clon,oid = case
    dargv = {} 
    dargv['oid'] = oid
    dargv['Year'] = Year
    dargv['Mon']  = Mon
    dargv['Day']  = Day

    dargv['iy']   = -9999
    dargv['ey']   = -9999
    dargv['clat'] = clat
    dargv['clon'] = clon
    dargv['DB_MAXREC'] = 20000
    dargv['dlatlon']   = 5
    dargv['xpos']      = 100

    if oid != 3556: continue # test

    #--------------------
    sargv = ['%s=%s'%(key, dargv[key]) for key in dargv.keys()]
    sargv = ' '.join(sargv)

    #-- Intersection (precipitaiton water content) ---
    '''
    prog = '/home/utsumi/bin/PMM/ret-epc/draw.intersect.prwat.py'
    lcmd = ['python', prog, sargv]
    print lcmd
    subprocess.call(lcmd)
    '''

    #-- Intersection (Zm) -----------------------------
    prog = '/home/utsumi/bin/PMM/ret-epc/draw.intersect.zm.py'
    lcmd = ['python', prog, sargv]
    print lcmd
    subprocess.call(lcmd)
 
 
