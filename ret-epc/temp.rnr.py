def read_rnr_table(idx_db):
    '''
    Rain/No Rain screening
    --- Contents of files ----
    calc 0:Skip all pixels  1:Screening  2:Calc all pixels
    D :  Linear discriminator threshold.
         If discriminator for the pixelsi smaller than this threshold, assign 'no-rain' to the pixel.
    SR:  Skip Ratio (expected fraction of skipped pixels)
    RMA: Ratio of Missing Amount  (>=0, 0.1, 1, 5, 10mm/h)
    m_org-s_no: First 12: For 12 EPCs. 13th: T2m
    '''

    rnrDir  = dbDir + '/rnr'
    srcPath = rnrDir + '/rnr.%05d.txt'%(idx_db)
    f=open(srcPath,'r'); lines=f.readlines(); f.close()

    rnrflag = int(lines[0].split('\t')[1])
    thD     = float(lines[1].split('\t')[1])
    SR      = float(lines[2].split('\t')[1])
    lRMA    = map(float, lines[3].split('\t')[1:])
    ave_org   = map(float, lines[4].split('\t')[1:])  # for normalization
    ave_rain  = map(float, lines[5].split('\t')[1:])  # Averages of normalized variables
    ave_no    = map(float, lines[6].split('\t')[1:])  # Averages of normalized variables
    std_org   = map(float, lines[7].split('\t')[1:])  # for normalization
    std_rain  = map(float, lines[8].split('\t')[1:])  # Averages of normalized variables
    std_no    = map(float, lines[9].split('\t')[1:])  # Averages of normalized variables

    return rnrflag, thD, SR, lRMA, ave_org, ave_rain, ave_no, std_org, std_rain, std_no

idx_db= 9865
dbDir = '/media/disk2/share/PMM/EPCDB/samp.20000.GMI.V05A.S1.ABp103-117.01-12'
l=read_rnr_table(idx_db)
for x in l:
    print x
