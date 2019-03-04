#! /usr/bin/python
#--------------------------------------------------------------------
# PROGRAM    : CC.Time-Height.py	
# CREATED BY : hjkim @IIS.2019-02-22 16:21:12.353039
# MODIFED BY :
#
# USAGE      : $ ./CC.Time-Height.py
#
# DESCRIPTION:
#------------------------------------------------------cf0.2@20120401


import  os,sys
from    optparse                import OptionParser

import  pylab   as pl
import  numpy   as np


def readdata( preset ):

    srcpath     = {'rmseNO': './data/table.nt-nlev.L2A25.5.0km.minNum.3.rmse.NO.RH.all.csv',
                   'rmseCR': './data/table.nt-nlev.L2A25.5.0km.minNum.3.rmse.CR.RH.all.csv',
                   'bratNO': './data/table.nt-nlev.L2A25.5.0km.minNum.3.brat.NO.RH.all.csv',
                   'bratCR': './data/table.nt-nlev.L2A25.5.0km.minNum.3.brat.CR.RH.all.csv',
                   'ccNO'  : './data/table.nt-nlev.L2A25.5.0km.minNum.3.cc.NO.RH.all.csv',
            }[ preset ]

    vmin, vmax, \
    cmap        = {'rmseNO': [ 3.8, 4.7,  pl.cm.autumn_r], 
                   'rmseCR': [ 3.8, 4.7,  pl.cm.autumn_r],
                   'bratNO': [-0.4, 0.4,  pl.cm.coolwarm],
                   'bratCR': [-0.4, 0.4,  pl.cm.coolwarm],
                   #'rmseNO': [ 3.7, 4.8,  pl.cm.Spectral_r], 
                   #'rmseCR': [ 3.7, 4.8,  pl.cm.Spectral_r],
                   #'bratNO': [-0.4, 0.1,  pl.cm.Spectral_r],
                   #'bratCR': [-0.4, 0.1,  pl.cm.Spectral_r],
                   'ccNO'  : [ 0.48,0.61, pl.cm.Spectral_r],
            }[ preset ]


    strdata     = np.array( [ l.strip().split(',') for l in open( srcpath ) ] )
    data        = strdata[1:,1:].astype('float32')

    xtlabels= strdata[0,1:]
    ytlabels= strdata[1:,0]

    normdata    = ( data - vmin ) / ( vmax - vmin )
    normdata    = (normdata * 100).astype('int8')

    print srcpath

    return data, normdata, xtlabels, ytlabels, vmin, vmax, cmap



def main(args,opts):
    print args
    print opts

    preset  = 'bratCR'
    preset  = 'rmseCR'
    preset  = 'bratNO'
    preset  = 'rmseNO'
    preset  = 'ccNO'

    nX      = 7
    nY      = 18

    X       = np.linspace( 0.08, 0.85, nX+1 )[:-1]
    Y       = np.linspace( 0.06, 0.94, nY+1 )[:-1]
    Y[0]   -= 0.02

    print X
    print Y

    width   = X[1]-X[0]-0.01
    height  = Y[2]-Y[1]-0.01

    fig     = pl.figure( figsize=(4.5,7) )
    ax0     = fig.add_axes([0,0,1,1])
    ax0.set_axis_off()
    ax0.text( 0.04, 0.5, 'Height', ha='center', va='center', rotation=90, fontsize=13 )
    ax0.text( 0.5, 0.02, 'Averaging Time', ha='center', va='center', fontsize=13 )

    Ax      = [ [ fig.add_axes( [ x,y,width,height ] ) for x in X ] 
            for y in Y ]

    data, normdata, xtlabels, ytlabels, vmin, vmax, cmap    = readdata( preset )
    print xtlabels
    print ytlabels
    print data.shape, data.min(), data.max(), vmin, vmax

    colors  = cmap( np.linspace(0,1, 101) )     

    print data.shape
    print len(Ax)

    for x, xtl in zip (X, xtlabels):
        ax0.text( x+0.05, 0.95, xtl, ha='center', va='center', fontsize=11 ) 

    for j, row in enumerate( Ax ):

        ax0.text( 0.86, Y[j]+0.02, ytlabels[j], ha='left', va='center', fontsize=11 )

        for i, ax in enumerate( row ):
            nd  = normdata[j,i]
            d   = data[j,i]
            ax.set_facecolor( colors[ nd ] )
            #ax.set_alpha(0.7)
            ax.get_xaxis().set_ticks([])
            ax.get_yaxis().set_ticks([])
            ax.text( 0.5,0.5, '{:3.2f}'.format( d ), ha='center', va='center', fontsize=11 )


    pl.savefig( 'Utsumi2019.{}.png'.format( preset ) )

    pl.show()

    return


if __name__=='__main__':
    usage   = 'usage: %prog [options] arg'
    version = '%prog 1.0'

    parser  = OptionParser(usage=usage,version=version)

#    parser.add_option('-r','--rescan',action='store_true',dest='rescan',
#                      help='rescan all directory to find missing file')

    (options,args)  = parser.parse_args()

#    if len(args) == 0:
#        parser.print_help()
#    else:
#        main(args,options)

#    LOG     = LOGGER()
    main(args,options)

    
