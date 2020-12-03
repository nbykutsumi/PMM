import glob

ssearch = '/home/utsumi/mnt/lab_tank/utsumi/PMM/MRMS/level2-pixel-match/*/*/*/30/num-all*'
print ssearch
lpath = sorted(glob.glob(ssearch))
for spath in lpath:
    print spath
