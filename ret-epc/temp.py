import numpy as np
from numpy import *

def read_table(srcPath, type=float):
    f=open(srcPath,'r'); lines=f.readlines(); f.close()
    lout = []
    for line in lines:
        line = line.strip().split()
        #line = map(type, line)
        lout.append(line)
    #return array(lout)
    return array(lout).astype(type)


coefDir  = '/media/disk2/share/PMM/JPLDB/EPC_COEF/GMI'
coefPath = coefDir + '/coef_pc.txt'
a2coef   = read_table(coefPath)
print (a2coef)
a2coef   = a2coef[:,1:]
#print a2coef
