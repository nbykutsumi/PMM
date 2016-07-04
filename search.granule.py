import gpm
import os
from datetime import datetime

prjName = 'GPM.KuPR.L2'
prdVer  = '02'

#prjName = 'GPM.KaPR.L2'
#prdVer  = '02'

#prjName = 'GPM.GMI.L2'
#prdVer  = '02'

BBox    = [[20.0, 118.0], [48.0, 150.0]]

iyear,eyear = 2014, 2014
imon, emon  = 4, 7
lyear       = range(iyear, eyear+1)
lmon        = range(imon,  emon+1)
for Y in lyear:
  for M in lmon:
    iY, iM = Y, M
    if M != 12:
      eY = iY
      eM = iM +1

    elif M == 12:
      eY = iY +1
      eM = 1 
    
    iDTime      = datetime(iY,iM,1,0)
    eDTime      = datetime(eY,eM,1,0)
    #eDTime      = datetime(eY,iM,1,5)
    
    OBgpm     = gpm.GPM(prjName, prdVer)
    lpath   = OBgpm.search_granules(BBox, iDTime, eDTime)
    
    sout  = "\n".join([s.split("/")[-1] for s in lpath]).strip()
    odir  = "/tank/utsumi/PMM/COMP.ORB.ORG/GranuleList/%s"%(prjName)
    
    try:
      os.makedirs(odir)
    except:
      pass
    
    oname = odir + "/Glist.LAT%05.2f-%05.2f.LON%05.2f-%05.2f.%04d.%02d.txt"%(BBox[0][0],BBox[1][0], BBox[0][1], BBox[1][1], Y,M)
    f=open(oname, "w"); f.write(sout); f.close()
    print oname
