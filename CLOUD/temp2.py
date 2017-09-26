from numpy import *
import myfunc.IO.GSMaP as GSMaP
from datetime import datetime 


DTime = datetime(2014,6,1,0)
gs = GSMaP.GSMaP(prj="standard",ver="v6", compressed=False)

a=gs.load_mmh(DTime)
print a
