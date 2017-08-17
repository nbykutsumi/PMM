import myfunc.util as util

Year = 1999
Mon  = 11
ldMon = range(30)

for dMon in ldMon:
    dMon = -dMon
    oYear,oMon = util.shift_YM(Year,Mon,dMon)
    print Year,Mon,"+",dMon,"-->",oYear,oMon
 
