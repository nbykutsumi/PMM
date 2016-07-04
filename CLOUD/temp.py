import myfunc.util as util
import myfunc.IO.GSMaP as GSMaP
from datetime import datetime

gsmap = GSMaP.GSMaP()
DTime = datetime(2015,4,3)
print gsmap.ret_path_sateinfo(DTime)
a   = gsmap.load_sateinfo(DTime)
print a

