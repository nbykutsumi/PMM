import GPMGV



gv = GPMGV.GPMGV()

#gv.load_sitelist_reclassified()



d = gv.ret_ddomYM2gName()
print d.keys()
key = d.keys()[0]
print ''
print key, d[key]
