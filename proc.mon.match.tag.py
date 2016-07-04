from numpy import *
from tag_func import *
import ctrack_func, tag_func, pmm_func
import ctrack_para, pmm_para
import ctrack_fig

calcflag = True
#calcflag = False
#figflag  = True
figflag  = False
iY, eY = 2002, 2009
#iY, eY = 2002, 2002
lY     = range(iY, eY+1)
lM = [1,2,3,4,5,6,7,8,9,10,11,12]

prtype1 = "GSMaP.v5"
prtype2 = "PR.2A25"
miss    = -9999.

para_tag = pmm_para.Para_tag()
ny,nx    = pmm_para.ret_nynx(prtype1)

para_tag.thpr = 0.0
para_tag.ny   = ny
para_tag.nx   = nx
para_tag.wnflag = False
if para_tag.wnflag == False:
  #----------
  ltag      = ["plain","tc","c","fbc","te","tf","ef","tef","ot"]
#  ltag      = ["plain"]
  ltag_comb = ["tc.comb","cf.comb"]

def ret_la2dat(Y,M,prtype1,prtype2, tag):
  sdir_root, sdir, sumname1, sumname2, numname \
      =tag_func.ret_name_tagthpr_match_sumnum_mon\
      (\
         prtype1, prtype2, tag, para_tag, Y, M\
      )
  a2sum1 = fromfile(sumname1, float32).reshape(ny,nx) 
  a2sum2 = fromfile(sumname2, float32).reshape(ny,nx) 
  a2num  = fromfile(numname,  float32).reshape(ny,nx) 
  return a2sum1, a2sum2, a2num

def mk_la2comb(tag_comb):
  if tag_comb == "tc.comb":
    a2sum_tc1,  a2sum_tc2 , a2num_tc  = ret_la2dat(Y,M,prtype1,prtype2,"tc")
    a2sum_te1,  a2sum_te2 , a2num_te  = ret_la2dat(Y,M,prtype1,prtype2,"te")
    a2sum_tf1,  a2sum_tf2 , a2num_tf  = ret_la2dat(Y,M,prtype1,prtype2,"tf")
    a2sum_tef1, a2sum_tef2, a2num_tef = ret_la2dat(Y,M,prtype1,prtype2,"tef")
    a2sum_pln1, a2sum_pln2, a2num_pln = ret_la2dat(Y,M,prtype1,prtype2,"plain")

    a2sum1 = a2sum_tc1 + a2sum_te1*0.5 + a2sum_tf1*0.5 + a2sum_tef1*(1./3.)
    a2sum2 = a2sum_tc2 + a2sum_te2*0.5 + a2sum_tf2*0.5 + a2sum_tef2*(1./3.)

    a2num_tag  = a2num_tc + a2num_te + a2num_tf + a2num_tef


    a2amt1  = (ma.masked_where(a2num_pln==0., a2sum1)/a2num_pln).filled(miss)
    a2amt2  = (ma.masked_where(a2num_pln==0., a2sum2)/a2num_pln).filled(miss)

    a2rat1  = (ma.masked_where(a2num_tag==0., a2sum1)/a2num_tag).filled(miss)
    a2rat2  = (ma.masked_where(a2num_tag==0., a2sum2)/a2num_tag).filled(miss)

  elif tag_comb == "cf.comb":
    a2sum_c1,  a2sum_c2  ,a2num_c   = ret_la2dat(Y,M,prtype1,prtype2,"c")
    a2sum_fbc1,a2sum_fbc2,a2num_fbc = ret_la2dat(Y,M,prtype1,prtype2,"fbc")
    a2sum_ef1, a2sum_ef2 ,a2num_ef  = ret_la2dat(Y,M,prtype1,prtype2,"ef")
    a2sum_te1, a2sum_te2 ,a2num_te  = ret_la2dat(Y,M,prtype1,prtype2,"te")
    a2sum_tf1, a2sum_tf2 ,a2num_tf  = ret_la2dat(Y,M,prtype1,prtype2,"tf")
    a2sum_tef1,a2sum_tef2,a2num_tef = ret_la2dat(Y,M,prtype1,prtype2,"tef")
    a2sum_pln1,a2sum_pln2,a2num_pln = ret_la2dat(Y,M,prtype1,prtype2,"plain")


    a2sum1 = a2sum_c1 + a2sum_fbc1 + a2sum_ef1\
            +a2sum_te1*0.5+ a2sum_tf1*0.5 + a2sum_tef1*(1./3.) 

    a2sum2 = a2sum_c2 + a2sum_fbc2 + a2sum_ef2\
            +a2sum_te2*0.5+ a2sum_tf2*0.5 + a2sum_tef2*(1./3.) 

    a2num_tag = a2num_c + a2num_fbc + a2num_ef\
               +a2num_te + a2num_tf + a2num_tef

    a2amt1  = (ma.masked_where(a2num_pln==0., a2sum1)/a2num_pln).filled(miss)
    a2amt2  = (ma.masked_where(a2num_pln==0., a2sum2)/a2num_pln).filled(miss)

    a2rat1  = (ma.masked_where(a2num_tag==0., a2sum1)/a2num_tag).filled(miss)
    a2rat2  = (ma.masked_where(a2num_tag==0., a2sum2)/a2num_tag).filled(miss)


  return a2amt1, a2amt2, a2rat1, a2rat2, a2sum1, a2sum2, a2num_tag

#***************************************************
if calcflag == True:
  for Y,M in [[Y,M] for Y in lY for M in lM]:
    #---- num plain ---------
    tag = "plain"
    a2num_pln = ret_la2dat(Y,M,prtype1,prtype2,"plain")[2]
    #---- calc precip (mean amount and rate) --
    for tag in ltag + ltag_comb:
      print ""
      print tag
      #------------
      if tag in ltag:
        a2sum1, a2sum2, a2num_tag\
         = ret_la2dat(Y,M,prtype1,prtype2,tag)

        a2amt1 = (ma.masked_where(a2num_pln==0., a2sum1)/a2num_pln).filled(miss)
        a2amt2 = (ma.masked_where(a2num_pln==0., a2sum2)/a2num_pln).filled(miss)
  
        a2rat1 = (ma.masked_where(a2num_tag==0., a2sum1)/a2num_tag).filled(miss)
        a2rat2 = (ma.masked_where(a2num_tag==0., a2sum2)/a2num_tag).filled(miss)

      elif tag in ltag_comb:
        a2amt1, a2amt2, a2rat1, a2rat2, a2sum1, a2sum2, a2num_tag = mk_la2comb(tag)

      #-- write ----
      sdir_root, sdir, sumpath1, sumpath2, numpath \
          =tag_func.ret_name_tagthpr_match_sumnum_mon\
          (\
             prtype1, prtype2, tag, para_tag, Y, M\
          )

      sname1 = sumpath1.split("/")[-1][4:]
      sname2 = sumpath2.split("/")[-1][4:]

      amtpath1 = sdir + "/" + "pr."+ sname1 
      amtpath2 = sdir + "/" + "pr."+ sname2 
      ratpath1 = sdir + "/" + "rat."+ sname1
      ratpath2 = sdir + "/" + "rat."+ sname2

      ctrack_func.mk_dir(sdir)
      a2sum1.tofile(sumpath1)
      a2sum2.tofile(sumpath2)
      a2num_tag.tofile(numpath)

      a2amt1.tofile(amtpath1)
      a2amt2.tofile(amtpath2)
      a2rat1.tofile(ratpath1)
      a2rat2.tofile(ratpath2)

      print amtpath1

 

