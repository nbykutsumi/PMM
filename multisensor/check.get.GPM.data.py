import glob

basedir = '/media/disk2/data/PMM/NASA'
ssearch = '/media/disk2/data/PMM/NASA/*/2A/V05/2018/*'

lpath = sorted(glob.glob(ssearch))
for spath in lpath:
    print spath
