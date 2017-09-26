from numpy import *
import myfunc.IO.IMERG as IMERG
from datetime import datetime 


DTime = datetime(2014,6,1,0)
imerg = IMERG.IMERG(PRD="PROD",VER="V04", crd="sa")

a=imerg.load_mmh(DTime)
print a
