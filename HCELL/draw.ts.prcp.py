from numpy import *
import os, sys
import myfunc.util          as util
import myfunc.IO.GPyML3     as GPyML3
import matplotlib.pyplot    as plt
import matplotlib

matplotlib.use('Agg')

iY  = 1999
eY  = 2002
#eY  = 2013
iYM =[iY,1]
eYM =[eY,12]
lYM = util.ret_lYM(iYM,eYM)
nYear = eY - iY +1

GRIDTYPE = 2
version  = "07"
figDir   = "/tank/utsumi/PMM/HCELL/pict"

gpm = GPyML3.L3A25(version=version, GRIDTYPE=GRIDTYPE, crd="sa")

ny  = gpm.ny
nx  = gpm.nx
Lat = gpm.Lat
Lon = gpm.Lon

lzlabel = ["2km","4km","6km","PathAve"]
# Functions ------------
def DrawTimeLat(a2in, figPath, stitle, vmin=None, vmax=None):
    fig = plt.figure(figsize=(6,3))
    ax  = fig.add_axes([0.1,0.15,0.85,0.64])
   
    ny,nx = a2in.shape 
    x   = range(nx)
    y   = Lat
    X,Y = meshgrid(x,y)

    im  = ax.pcolormesh(X,Y,a2in,vmin=vmin,vmax=vmax)

    # Colorbar
    cb  = plt.colorbar(im, orientation="horizontal")
    cb.ax.tick_params(labelsize=9)

    # X-tick label
    plt.xticks(x[::12])
    lxlabel = range(iYM[0],eYM[0]+1)
    ax.xaxis.set_ticklabels(lxlabel, fontsize=8)

    # Y-tick label
    yticks  = arange(-35,35+0.1,5)
    plt.yticks(yticks)
    ax.yaxis.set_ticklabels(yticks, fontsize=8)

    # Title
    plt.title(stitle,fontsize=9)

    # Save
    plt.savefig(figPath)
    print figPath


def mkClimMon(a2in):
    # Seasonal climatology
    a2clim    = empty(a2in.shape)
    for iM in range(12):
        a1climMon= a2in[:,iM::12].mean(axis=1).reshape(ny,-1)
        a2clim[:,iM::12] = a1climMon
    return a2clim


#-----------------------
"""
#*****************************
# Non-zero rain rate
Var = "rainMean2"
nz  = 4
da2dat = {i:empty([ny,len(lYM)]) for i in range(nz)}
for itime, (Y,M) in enumerate(lYM):
    a3in = gpm.load_var(Y,M,Var)
    for iz in range(nz):
        da2dat[iz][:,itime] = ma.masked_less(a3in[iz],0).mean(axis=1)


# Figure
for iz in range(nz):
    figPath = os.path.join(figDir, "TimeLat.%s.%02d.png"%(Var, iz))


    stitle  = "non-zero rain rate 3A25 %s @%s [mm/h]"\
                %(Var, lzlabel[iz])
    DrawTimeLat(da2dat[iz], figPath, stitle)
"""


"""
#*****************************
# Non-zero rain rate (Normarized)
Var = "rainMean2"
nz  = 4
da2dat = {i:empty([ny,len(lYM)]) for i in range(nz)}
da2out = {i:empty([ny,len(lYM)]) for i in range(nz)}
for itime, (Y,M) in enumerate(lYM):
    a3in = gpm.load_var(Y,M,Var)
    for iz in range(nz):
        da2dat[iz][:,itime] = ma.masked_less(a3in[iz],0).mean(axis=1)

# Normalize
for iz in range(nz):
    a2clim     = mkClimMon(da2dat[iz])
    da2out[iz] = (da2dat[iz] - a2clim)/a2clim


# Figure
for iz in range(nz):
    figPath = os.path.join(figDir, "TimeLat.%s.Norm.%02d.png"%(Var, iz))


    stitle  = "normalized non-zero rain rate 3A25 %s @%s [mm/h]"\
                %(Var, lzlabel[iz])
    DrawTimeLat(da2out[iz], figPath, stitle)
"""

"""
#*****************************
# Conv rain rate
nz  = 1
Var = "surfRainConvMean2"
da2dat = {i:empty([ny,len(lYM)]) for i in range(nz)}
for itime, (Y,M) in enumerate(lYM):
    a3in = gpm.load_var(Y,M,Var)
    for iz in range(nz):
        da2dat[iz][:,itime] = ma.masked_less(a3in[iz],0).mean(axis=1)


# Figure
for iz in range(nz):
    figPath = os.path.join(figDir, "TimeLat.%s.%02d.png"%(Var, iz))


    stitle  = "conv near surf rain 3A25 %s [mm/h]"\
                %(Var)
    DrawTimeLat(da2dat[iz], figPath, stitle)

"""


"""
#*****************************
# Conv rain rate (Normalized)
Var = "surfRainConvMean2"
nz  = 1
da2dat = {i:empty([ny,len(lYM)]) for i in range(nz)}
da2out = {i:empty([ny,len(lYM)]) for i in range(nz)}
for itime, (Y,M) in enumerate(lYM):
    a3in = gpm.load_var(Y,M,Var)
    for iz in range(nz):
        da2dat[iz][:,itime] = ma.masked_less(a3in[iz],0).mean(axis=1)

# Normalize
for iz in range(nz):
    a2clim     = mkClimMon(da2dat[iz])
    da2out[iz] = (da2dat[iz] - a2clim)/a2clim


# Figure
for iz in range(nz):
    figPath = os.path.join(figDir, "TimeLat.%s.Norm.%02d.png"%(Var, iz))


    stitle  = "normalized near surf conv rain 3A25 %s [mm/h]"\
                %(Var)
    DrawTimeLat(da2out[iz], figPath, stitle)
"""
"""
#*****************************
# Conv rain pixels
nz  = 1
Var = "surfRainConvPix2"
da2dat = {i:empty([ny,len(lYM)]) for i in range(nz)}
for itime, (Y,M) in enumerate(lYM):
    a3in = gpm.load_var(Y,M,Var)
    for iz in range(nz):
        da2dat[iz][:,itime] = ma.masked_less(a3in[iz],0).mean(axis=1)


# Figure
for iz in range(nz):
    figPath = os.path.join(figDir, "TimeLat.%s.%02d.png"%(Var, iz))


    stitle  = "conv near surf rain counts 3A25 %s [counts]"\
                %(Var)
    DrawTimeLat(da2dat[iz], figPath, stitle)
"""

"""
#*****************************
# Conv rain pixels (Normalized)
Var = "surfRainConvPix2"
nz  = 1
da2dat = {i:empty([ny,len(lYM)]) for i in range(nz)}
da2out = {i:empty([ny,len(lYM)]) for i in range(nz)}
for itime, (Y,M) in enumerate(lYM):
    a3in = gpm.load_var(Y,M,Var)
    for iz in range(nz):
        da2dat[iz][:,itime] = ma.masked_less(a3in[iz],0).mean(axis=1)

# Normalize
for iz in range(nz):
    a2clim     = mkClimMon(da2dat[iz])
    da2out[iz] = (da2dat[iz] - a2clim)/a2clim


# Figure
for iz in range(nz):
    figPath = os.path.join(figDir, "TimeLat.%s.Norm.%02d.png"%(Var, iz))


    stitle  = "normalized near surf conv counts 3A25 %s [counts]"\
                %(Var)
    DrawTimeLat(da2out[iz], figPath, stitle)
"""


"""
#*****************************
# Shallow rain pixels
nz  = 1
Var = "shallowRainPix2"
da2dat = {i:empty([ny,len(lYM)]) for i in range(nz)}
for itime, (Y,M) in enumerate(lYM):
    a3in = gpm.load_var(Y,M,Var)
    for iz in range(nz):
        da2dat[iz][:,itime] = ma.masked_less(a3in[iz],0).mean(axis=1)


# Figure
for iz in range(nz):
    figPath = os.path.join(figDir, "TimeLat.%s.%02d.png"%(Var, iz))


    stitle  = "shallow rain counts 3A25 %s [counts]"\
                %(Var)
    DrawTimeLat(da2dat[iz], figPath, stitle)
"""

"""
#*****************************
# Shallow rain pixels (Normalized)
Var = "shallowRainPix2"
nz  = 1
da2dat = {i:empty([ny,len(lYM)]) for i in range(nz)}
da2out = {i:empty([ny,len(lYM)]) for i in range(nz)}
for itime, (Y,M) in enumerate(lYM):
    a3in = gpm.load_var(Y,M,Var)
    for iz in range(nz):
        da2dat[iz][:,itime] = ma.masked_less(a3in[iz],0).mean(axis=1)

# Normalize
for iz in range(nz):
    a2clim     = mkClimMon(da2dat[iz])
    da2out[iz] = (da2dat[iz] - a2clim)/a2clim


# Figure
for iz in range(nz):
    figPath = os.path.join(figDir, "TimeLat.%s.Norm.%02d.png"%(Var, iz))


    stitle  = "normalized shallow counts 3A25 %s [counts]"\
                %(Var)
    DrawTimeLat(da2out[iz], figPath, stitle)
"""


"""
#*****************************
# Iso Shallow rain pixels
nz  = 1
Var = "shallowIsoRainPix2"
da2dat = {i:empty([ny,len(lYM)]) for i in range(nz)}
for itime, (Y,M) in enumerate(lYM):
    a3in = gpm.load_var(Y,M,Var)
    for iz in range(nz):
        da2dat[iz][:,itime] = ma.masked_less(a3in[iz],0).mean(axis=1)


# Figure
for iz in range(nz):
    figPath = os.path.join(figDir, "TimeLat.%s.%02d.png"%(Var, iz))


    stitle  = "Iso shallow rain counts 3A25 %s [counts]"\
                %(Var)
    DrawTimeLat(da2dat[iz], figPath, stitle)
"""

"""
#*****************************
# Iso Shallow rain pixels (Normalized)
Var = "shallowIsoRainPix2"
nz  = 1
da2dat = {i:empty([ny,len(lYM)]) for i in range(nz)}
da2out = {i:empty([ny,len(lYM)]) for i in range(nz)}
for itime, (Y,M) in enumerate(lYM):
    a3in = gpm.load_var(Y,M,Var)
    for iz in range(nz):
        da2dat[iz][:,itime] = ma.masked_less(a3in[iz],0).mean(axis=1)

# Normalize
for iz in range(nz):
    a2clim     = mkClimMon(da2dat[iz])
    da2out[iz] = (da2dat[iz] - a2clim)/a2clim


# Figure
for iz in range(nz):
    figPath = os.path.join(figDir, "TimeLat.%s.Norm.%02d.png"%(Var, iz))


    stitle  = "normalized Iso shallow counts 3A25 %s [counts]"\
                %(Var)
    DrawTimeLat(da2out[iz], figPath, stitle)
"""

#"""
#*****************************
# storm height 
nz  = 2
Var = "stormHeightMean"
lzname= ["strat","conv"]
da2dat = {i:empty([ny,len(lYM)]) for i in range(nz)}
for itime, (Y,M) in enumerate(lYM):
    a3in = gpm.load_var(Y,M,Var)
    for iz in range(nz):
        da2dat[iz][:,itime] = ma.masked_less(a3in[iz],0).mean(axis=1)


# Figure
for iz in range(nz):
    figPath = os.path.join(figDir, "TimeLat.%s.%02d.png"%(Var, iz))


    stitle  = "Storm Height 3A25 %s [m] @%s"\
                %(Var, lzname[iz])
    vmin = 0
    vmax = 6000
    DrawTimeLat(da2dat[iz], figPath, stitle, vmin, vmax)
#"""

"""
#*****************************
# storm height (Normalized)
nz  = 2
Var = "stormHeightMean"
lzname= ["strat","conv"]
da2dat = {i:empty([ny,len(lYM)]) for i in range(nz)}
da2out = {i:empty([ny,len(lYM)]) for i in range(nz)}
for itime, (Y,M) in enumerate(lYM):
    a3in = gpm.load_var(Y,M,Var)
    for iz in range(nz):
        da2dat[iz][:,itime] = ma.masked_less(a3in[iz],0).mean(axis=1)

# Normalize
for iz in range(nz):
    a2clim     = mkClimMon(da2dat[iz])
    da2out[iz] = (da2dat[iz] - a2clim)/a2clim


# Figure
for iz in range(nz):
    figPath = os.path.join(figDir, "TimeLat.%s.Norm.%02d.png"%(Var, iz))


    stitle  = "normalized Storm Height 3A25 %s [m] @%s"\
                %(Var,lzname[iz])
    DrawTimeLat(da2out[iz], figPath, stitle)
"""


