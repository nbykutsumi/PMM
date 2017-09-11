from  numpy import *
import sys

def ret_regionBBox(region):
    # Southern Hemisphere
    BBox = None
    lllat = 0
    urlat = 37
    if region=="NAT":
        BBox = [[lllat,-75],[urlat,-15]]
    if region=="NAF":
        BBox = [[lllat,-15],[urlat,50]]
    if region=="ASI":
        BBox = [[lllat,50],[urlat,125]]
    if region=="NPA":
        BBox = [[lllat,125],[urlat,-125]]
    if region=="NAM":
        BBox = [[lllat,-125],[urlat,-75]]

    lllat = -37
    urlat = 0
    if region=="SAT":
        BBox = [[lllat,-35],[urlat,10]]
    if region=="SAF":
        BBox = [[lllat,10],[urlat,45]]
    if region=="SIN":
        BBox = [[lllat,45],[urlat,110]]
    if region=="OCE":
        BBox = [[lllat,110],[urlat,155]]
    if region=="SWP":
        BBox = [[lllat,155],[urlat,-150]]
    if region=="SEP":
        BBox = [[lllat,-150],[urlat,-80]]
    if region=="SAM":
        BBox = [[lllat,-80],[urlat,-35]]

    if BBox==None:
        print "check BBox",BBox
        sys.exit()
    return BBox


#def ret_regionMask(code):
#        a2region= grids.mk_mask_BBox(Lat, Lon, BBox)
 
