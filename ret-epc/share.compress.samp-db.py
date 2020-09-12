import subprocess, os, sys

ibasedir = '/media/disk2/share/PMM/EPCDB/samp.10000.GMI.V05A.S1.ABp103-117.01-12'
obasedir = '/home/utsumi/mnt/lab_home_rainbow/public_html/outgoing/ret-epc-202003/EPCDB'
#obasedir = '/home/utsumi/work_space'
lvname = ['DPRGMI_NS_surfPrecipTotRate','DPRGMI_NS_vfracConv','Ka_MS_precipRateNearSurface']
#lvname = ['DPRGMI_NS_surfPrecipTotRate']

for vname in lvname:
    os.chdir(ibasedir)

    srcpath = './%s'%(vname)
    print os.path.exists(srcpath)
    outpath = obasedir + '/%s.tar.gz'%(vname)
    print outpath
    scmd = 'tar -zcvf %s %s'%(outpath, srcpath)
    lcmd = scmd.split()
    subprocess.call(lcmd)
    print scmd

