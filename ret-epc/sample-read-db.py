import JPLDB


db = JPLDB.JPLDB()
srcPath = '/work/hk01/utsumi/JPLDB/EPC_DB/GMI_EPC_DATABASE/db_14603.bin'

db.set_file(srcPath) # set file

#-- Read Tb --
Tb  = db.get_var('tb')

#-- Read pixel Lat and Lon --
Lat = db.get_var('glat')
Lon = db.get_var('glon')


print 'shape of Tb'
print Tb.shape
print 'shape of Lat'
print Lat.shape
print 'shape of Lon'
print Lon.shape


