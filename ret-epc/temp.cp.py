import shutil, os
import subprocess
import myfunc.util as util

oid = 23936

lvname = ['DPRGMI_NS_surfPrecipTotRate','Ka_MS_precipRateNearSurface','DPRGMI_NS_precipTotWaterContRelSurf','DPRGMI_NS_vfracConv']

basedir='/media/disk2/share/PMM/EPCDB/samp.10000.GMI.V05A.S1.ABp103-117.01-12'

for vname in lvname:
    srcpath = basedir + '/%s/%05d.npy'%(vname, oid)
    print srcpath
    print os.path.exists(srcpath) 

    odir  = '/home/utsumi/mnt/lab_home_rainbow/public_html/outgoing/ret-epc-202003/temp/%s'%(vname)
    util.mk_dir(odir)
    
    opath = odir + '/' + os.path.basename(srcpath)
    print opath

    shutil.copy(srcpath, opath)
