import os, sys
import myfunc.util as util
import socket

myhost = socket.gethostname()
if myhost =='shui':
    pairbaseDir = '/tank/utsumi/env/pair'
    figDir = '/home.rainbow/utsumi/public_html/tempfig'
elif myhost == 'well':
    pairbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/env/pair'
    figDir = '/home/utsumi/mnt/lab_home_rainbow/public_html/tempfig'
else:
    print 'check myhost'
    sys.exit()


lsurftype = ['sea','land']

lvar = ['w_1.5','w_4.5','w_7.5']
lvar = lvar+['r_1.5','r_4.5','r_7.5']
lvar = lvar+['mvimd']
lvar = lvar+['tcwv']
lvar = lvar+['cape']
lvar = lvar+['skt']
lvar = lvar+['deptdz_low','deptdz_mid']
lvar = lvar+['dtvdz_low','dtvdz_mid']
#lvar = ['pre.deptdz_low']
#lvar = ['deptdz_low','deptdz_mid']
#lvar = ['dtvdz_low','dtvdz_mid']

#lvar = ['w_1.5','w_4.5','w_7.5']
#lvar = lvar+['r_1.5','r_4.5','r_7.5']
#lvar = lvar+['mvimd','tcwv']
#lvar = lvar+['cape']
#lvar = lvar+['dtvdz_low','dtvdz_mid']
#lvar = lvar+['pre.deptdz_low','pre.deptdz_mid']
#lvar = lvar+['pre.dtvdz_low','pre.dtvdz_mid']
#lvar = lvar+['pre.skt']

#lvar = ['deptdz_low','deptdz_mid']
#lvar = ['pre.dtvdz_low','pre.dtvdz_mid']
#lvar = ['pre.cape']

for surftype in lsurftype:
    sout = ''
    for var in lvar:
        corrDir = figDir + '/txt'
        util.mk_dir(corrDir)
        corrPath= corrDir + '/cc.%s.%s.txt'%(var,surftype)
        f=open(corrPath,'r');cc=f.read().strip(); f.close()
        print surftype, var, cc
        sout = sout+ '%s,%s,%s\n'%(surftype,var,cc)

    txtPath = corrDir + '/cc.summary.%s.csv'%(surftype)
    f=open(txtPath, 'w'); f.write(sout); f.close()
    print txtPath

