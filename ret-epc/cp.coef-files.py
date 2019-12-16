import os, shutil

lsensor = ['AMSR2','SSMIS','MHS','ATMS']
#fileName = 'coef_pc.txt'
fileName = 'ave_pc.txt'
for sensor in lsensor:
    iDir = '/media/disk2/share/PMM/EPCDB/EPC_COEF/%s'%(sensor)
    iPath= iDir + '/%s'%(fileName)
    oDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/JPLDB/EPC_COEF/%s'%(sensor)

    print sensor, os.path.exists(iPath), os.path.exists(oDir)
    shutil.copy(iPath, oDir)

