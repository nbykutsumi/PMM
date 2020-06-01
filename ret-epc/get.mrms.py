# %%
from datetime import datetime
import myfunc.util as util
import calendar
import os, sys
import socket
from datetime import datetime, timedelta
import numpy as np
import subprocess

hostname  = "https://pmm-gv.gsfc.nasa.gov/"
myhost = socket.gethostname()
if myhost =="shui":
    orootDir   = "/work/hk02/PMM/MRMS/level2"
elif myhost =="well":
    orootDir   = "/home/utsumi/mnt/lab_work/hk02/PMM/MRMS/level2"

iYM = [2018,7]
eYM = [2018,12]
lYM = util.ret_lYM(iYM,eYM)

gmi       = ["GPM","GMI"]
amsr2     = ["GCOMW1","AMSR2"]
ssmis_f16 = ["F16","SSMIS"]
ssmis_f17 = ["F17","SSMIS"]
ssmis_f18 = ["F18","SSMIS"]
atms_npp  = ["NPP","ATMS"]
atms_noaa20= ["NOAA20","ATMS"]

mhs_metopa= ["METOPA","MHS"]
mhs_metopb= ["METOPB","MHS"]
mhs_noaa18= ["NOAA18","MHS"]
mhs_noaa19= ["NOAA19","MHS"]

lspec = [gmi, amsr2, ssmis_f16,ssmis_f17,ssmis_f18,atms_npp, atms_noaa20, mhs_metopa, mhs_metopb, mhs_noaa18, mhs_noaa19]
#lspec = [amsr2, ssmis_f16,ssmis_f17,ssmis_f18,atms_npp, atms_noaa20, mhs_metopa, mhs_metopb, mhs_noaa18, mhs_noaa19]

lvar = ['PRECIPRATE','MASK','RQI']
#----------------------------------
def mk_dir(sodir):
  try:
    os.makedirs(sodir)
  except OSError:
    pass


#----------------------------------
for (Year,Mon) in lYM:
    for spec in lspec:
        sate      = spec[0]
        sensor    = spec[1]

        if (sensor=='GMI') &(Mon<=7):continue
        if (sate=='GCOMW1')&(Mon<=7):continue
        if (sate=='F16')   &(Mon<=7):continue
        #if (sate=='F17')   &(Mon<=4):continue
        #if (sate=='F18')   &(Mon<=4):continue
        #if (sate=='NPP')   &(Mon<=3):continue
        #if (sate=='NOAA20')&(Mon<=3):continue
        #if (sate=='METOPA')&(Mon<=3):continue
        #if (sate=='METOPB')&(Mon<=3):continue
        #if (sate=='NOAA18')&(Mon<=3):continue
        #if (sate=='NOAA19')&(Mon<=3):continue

        #--- path and directory: Remote -----------
        irootDir = '/pub/NMQ/level2'
        iDir = irootDir + "/%s/%04d/%02d"%(sate,Year,Mon)
        oDir = orootDir + "/%s/%04d/%02d"%(sate,Year,Mon) 
        mk_dir(oDir)
        #--- list --------------

        for var in lvar:
            #if (sate=='METOPB')&(Mon==3)&(var in ['PRECIPRATE','MASK']):continue
            #iPath = iDir + '/' + '%s.*.gz'
            scmd  = 'wget -r --no-directories --no-parent -A %s.*.gz %s/%s/ --directory-prefix=%s'%(var, hostname, iDir, oDir) 
            #scmd  = 'wget -A  %s/%s/%s.*.gz --directory-prefix=%s'%(hostname,iDir,var,oDir) 
            lcmd  = scmd.split()

            subprocess.call(lcmd)





# %%
