import os, sys
import myfunc.util as util
import socket

myhost = socket.gethostname()
if myhost =='shui':
    pairbaseDir = '/tank/utsumi/env/pair'
    figDir = '/home.rainbow/utsumi/public_html/tempfig'
    figDir = '/home/utsumi/temp/env'
elif myhost == 'well':
    pairbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/env/pair'
    figDir = '/home/utsumi/mnt/lab_home_rainbow/public_html/tempfig'
    figDir = '/home/utsumi/temp/env'
else:
    print 'check myhost'
    sys.exit()


lsurftype = ['sea','land']
lheight   = ['all','low','high']
lts       = ['cold','cool','warm','hot']
lvar = ['r_15','r_45','r_75']
lvar = lvar+['cape']
lvar = lvar+['skt']
lvar = lvar+['dept_low','dept_mid']
lvar = lvar+['dtv_low','dtv_mid']
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
season = 7

sout = ''
for surftype in lsurftype:
    for height in lheight:
        for ts in lts:
            for var in lvar:
                corrDir = figDir + '/txt'
                util.mk_dir(corrDir)
                corrPath= corrDir + '/cc.conv.%s.dt-01.%s.%s.%s.%s.txt'%(var,season,height,ts,surftype)
                f=open(corrPath,'r');cc=f.read().strip(); f.close()
                print surftype, var, cc
                sout = sout+ '%s,%s,%s,%s,%s\n'%(surftype,height, ts, var, cc)

txtPath = corrDir + '/cc.summary.csv'
f=open(txtPath, 'w'); f.write(sout); f.close()
print txtPath

