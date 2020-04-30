import os, shutil
import glob

lvname = ['qualityFlag','precipTotWaterCont']
lmon = range(1,12+1)

basedir0 = '/media/disk2/share/PMM/EPCDB/samp.10000.GMI.V05A.S1.ABp103-117.01-12'

basedir1 = '/home/utsumi/mnt/lab_tank/utsumi/PMM/EPCDB/samp.10000.GMI.V05A.S1.ABp103-117.01-12'


for basedir in [basedir0, basedir1]:
    lvdir = glob.glob(basedir + '/*')

    for vdir in lvdir:
        for mon in lmon:
            srcdir = vdir + '/2017%02d'%(mon)
            #print os.path.exists(srcdir)
            if os.path.exists(srcdir):
                print srcdir

#for vname in lvname:
#    for basedir in [basedir0, basedir1]:
#        for mon in lmon:
#            srcdir = basedir + '/%s/2017%02d'%(vname,mon)
#            print ''
#            print vname, os.path.exists(srcdir)
#            print srcdir
#            #if os.path.exists(srcdir):
#            #    shutil.rmtree(srcdir)
