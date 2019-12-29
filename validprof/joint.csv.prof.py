import numpy as np
import csv
import pandas as pd

lprrange=[[0.5,999]]
lstype = ['veg','sea','snow']
lptype = ['conv','stra']
lseason=['ALL','DJF','JJA']
lregion = ['MIDN','TRO']
srcDir = '/home/utsumi/temp/ret'
lph    = ['A']
thpr0,thpr1 = 0.5, 999
for ph in lph:
    l = []
    for region in lregion:
        for season in lseason:
            if (region=='MIDN')and(season=='ALL'): continue
            if (region=='TRO')and(season in ['DJF','JJA']): continue
    
            for stype in lstype:
                for ptype in lptype:
                    
                    stampOut  = 's-%s.p-%s.ph-%s.pr-%.1f-%.1f.%s'%(stype,ptype,ph,thpr0,thpr1,season)
                    csvPath = srcDir + '/prof.%s.%s.csv'%(stampOut,region)
                    l.append(pd.read_csv(csvPath))
    
    df = pd.concat(l, axis=1)
    outPath = srcDir + '/joint.prof.%s.csv'%(ph)
    df.to_csv(outPath)
    print outPath
