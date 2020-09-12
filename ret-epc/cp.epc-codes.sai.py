import os, sys, shutil
import myfunc.util as util
import subprocess

lfname=[
'/home/utsumi/bin/PMM/ret-epc/ret-epc-multi.2020.07.20.py',
'/home/utsumi/bin/PMM/ret-epc/run-epc-simple.2020.08.17.py',
'/home/utsumi/bin/PMM/ret-epc/JPLDB.py',
'/home/utsumi/bin/PMM/ret-epc/EPCDB.py',
'/home/utsumi/bin/PMM/ret-epc/epcfunc.py',
'/home/utsumi/bin/PMM/ret-epc/README',
'/home/utsumi/mnt/lab_work/hk02/PMM/NASA/GPM.GMI/1C/V05/2014/06/01/1C.GPM.GMI.XCAL2016-C.20140601-S010534-E023807.001453.V05A.HDF5',
'/home/utsumi/mnt/lab_tank/utsumi/PMM/MATCH.GMI.GPM.V05/S1.ABp000-220.GMI.S2.IDX/2014/06/01/Xpy.1.001453.npy',
'/home/utsumi/mnt/lab_tank/utsumi/PMM/MATCH.GMI.GPM.V05/S1.ABp000-220.GMI.S2.IDX/2014/06/01/Ypy.1.001453.npy',
'/home/utsumi/bin/PMM/mkdb/f_match_fov.f90',
'/home/utsumi/bin/PMM/mkdb/f2py.make.py',
'/home/utsumi/bin/PMM/mkdb/mk.match.idx.gmiS2.gmi.fullswath.py',
'/home/utsumi/mnt/lab_tank/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.GMI.S2.dydx/dx.000.npy',
'/home/utsumi/mnt/lab_tank/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.GMI.S2.dydx/dy.000.npy',
'/home/utsumi/mnt/lab_tank/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.GMI.S2.dydx/dx.180.npy',
'/home/utsumi/mnt/lab_tank/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.GMI.S2.dydx/dy.180.npy',
'/home/utsumi/mnt/lab_tank/utsumi/PMM/MATCH.GMI.GPM.V05/S1.ABp000-220.gtopo/2014/06/01/gtopo.001453.npy',
'/home/utsumi/mnt/lab_tank/utsumi/PMM/MATCH.GMI.GPM.V05/S1.ABp000-220.MERRA2.t2m/2014/06/01/t2m.001453.npy',
'/home/utsumi/mnt/lab_tank/utsumi/PMM/MATCH.GMI.GPM.V05/S1.ABp000-220.MERRA2.tqv/2014/06/01/tqv.001453.npy',
]

distdir = '/home/utsumi/temp/share-epc/epc-2020-08-17'
util.mk_dir(distdir)

for fname in lfname:
    oname = distdir + '/' + os.path.basename(fname)
    shutil.copy(fname, oname)
    print oname

#-- make tar.gz ----
distrootdir = os.path.dirname(distdir)
tarname = distdir + '.tar.gz'
cmd = 'tar -zcvf %s -C %s %s'%(tarname, distrootdir, os.path.basename(distdir))
subprocess.call(cmd.split())
print cmd

#-- copy tar.gz to public_html----
outgoingdir = '/home/utsumi/mnt/lab_home_rainbow/public_html/outgoing/ret-epc-202003'
taroutname  = outgoingdir + '/' + os.path.basename(tarname)
shutil.copy(tarname, taroutname)
print taroutname
