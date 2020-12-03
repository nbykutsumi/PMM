import glob

ssearch = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retepc/us.v01/*/2018/*'
print ssearch
lpath = sorted(glob.glob(ssearch))
for spath in lpath:
    print spath
