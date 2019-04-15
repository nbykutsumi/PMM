import subprocess
import sys
import socket
hostname = socket.gethostname()

srcName = sys.argv[1]
if srcName[-2:]==".c":
    binName = srcName[:-2]
else:
    print "check file name",srcName
#-- compile -------
if hostname=='mizu':
    cmd = "gcc -fpack-struct %s -L/usr/local/lib -lnetcdf -lm -I/usr/local/include -o %s"%(srcName, binName)

elif hostname=='shui':
    #cmd = "gcc -fpack-struct %s -L/usr/local/lib -lnetcdf -lm -I/usr/local/include -o %s"%(srcName, binName)
    cmd = "gcc -fpack-struct %s -L/usr/local/lib -lnetcdf -lm -I/usr/local/include amsr2.o atms.o gmi.o mhs.o saphir.o ssmis.o -o %s"%(srcName, binName)

elif hostname=='well':
    cmd = "gcc -fpack-struct %s -L/usr/local/lib -lnetcdf -lm -I/usr/local/include -o %s"%(srcName, binName)

else:
    print 'check hostname'
    sys.exit() 

print cmd
subprocess.call(cmd, shell=True)
print "Compiled"
print binName

