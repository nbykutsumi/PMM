import os, sys, re
from numpy import *
import pmm_func
import calendar

lY = [2001]
lM = [1]
lH = range(0,23+1)

for Y,M in [[Y,M] for Y in lY for M in lM]:
  eD = calendar.monthrange(Y,M)[1]
  for D in range(1,eD+1):
    for H in lH:
      a2pr = pmm_func.ret_a2backward_org(Y,M,D,H,minute=0,prtype="GSMaP.v5",maskmiss=True)
      print Y,M,D,H,"**", a2pr.min(), a2pr.max(), a2pr.mean()
