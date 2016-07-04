from numpy import *
from tag_func import *
import ctrack_func, tag_func, pmm_func
import ctrack_para, pmm_para
import ctrack_fig

#lseason = ["DJF","JJA"]
#lseason = ["ALL"]
#lseason = ["MAM","DJF","JJA","SON","ALL"]
lseason = [1,2,3,4,5,6,7,8,9,10,11,12]
#lseason = [1]
iY, eY = 2002, 2009
#iY, eY = 2002, 2002
lY     = range(iY, eY+1)

prtype1 = "GSMaP.v5"
prtype2 = "PR.2A25"
miss    = -9999.
lvar    = ["pr","rat"]
#lvar    = ["rat"]
para_tag = pmm_para.Para_tag()
ny,nx    = pmm_para.ret_nynx(prtype1)

para_tag.thpr = 0.0
para_tag.ny   = ny
para_tag.nx   = nx
para_tag.wnflag = False
if para_tag.wnflag == False:
  #----------
  #ltag      = ["plain","tc","c","fbc","te","tf","ef","tef","ot","tc.comb","cf.comb"]
  ltag      = ["plain","ot","tc.comb","cf.comb"]

def ret_la2dat(var, Y,M,prtype1,prtype2, tag):
  sdir_root, sdir, varpath1, varpath2\
      =tag_func.ret_name_tagthpr_match_mon\
      (\
         var, prtype1, prtype2, tag, para_tag, Y, M\
      )
  a2var1 = fromfile(varpath1,  float32).reshape(ny,nx) 
  a2var2 = fromfile(varpath2,  float32).reshape(ny,nx) 
  return a2var1, a2var2



def ret_name_clim(prtype1,prtype2,tag,para_tag,iY,eY,season):
  lname = tag_func.ret_name_tagthpr_match_sumnum_clim\
          (\
           prtype1, prtype2, tag, para_tag, iY, eY, season\
          )

  sdir_root, sdir, amtname1, amtname2, ratname1, ratname2, sumname1, sumname2, numname\
        = lname
  return sdir_root, sdir, amtname1, amtname2, ratname1, ratname2, sumname1, sumname2, numname

def ret_a2sig(a2svar, a2ssqrt, nstep):
  N = nstep
  a2ave = a2svar / N
  a2sig = sqrt\
          (\
           (N * square(a2ave) - 2.0*a2ave * a2svar + a2svv)/N\
          ) 
  return a2sig 
      

#***************************************************
for var in lvar:
  for season in lseason:
    lM = ctrack_para.ret_lmon(season)
    #---- calc precip (mean amount and rate) --
    for tag in ltag:
    #for tag in ["cf.comb"]:  # test
      a2svar  = zeros([ny,nx],float32)
      a2svv   = zeros([ny,nx],float32)
  
      nstep  = 0

      #ltmp  = []
      for Y,M in [[Y,M] for Y in lY for M in lM]:
        nstep = nstep + 1
  
        a2var1, a2var2 = ret_la2dat(var,Y,M,prtype1,prtype2,tag)
        a2var1 = ma.masked_equal(a2var1, miss).filled(0.0)
        a2var2 = ma.masked_equal(a2var2, miss).filled(0.0)
        a2dvar = a2var1 - a2var2
  
        a2svar = a2svar + a2dvar
        a2svv  = a2svv  + square(a2dvar)

      #  ltmp.append(a2dvar[1100,1400])
      #----------
      a2sig = ret_a2sig(a2svar, a2svv, nstep)
      #-- write --
      sdir_root, sdir, oname, tmp = ret_name_tagthpr_match_clim("sig.d%s"%(var), prtype1, prtype2, tag, para_tag, iY, eY, season) 
  
      a2sig.tofile(oname)
      #print oname 
      #print "**************"
      #print "ltmp",ltmp
      #print var,std(ltmp), a2sig[1100,1400], mean(ltmp)

