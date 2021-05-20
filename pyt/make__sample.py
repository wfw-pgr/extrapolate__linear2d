import numpy as np


# ========================================================= #
# ===  make__sample.py                                  === #
# ========================================================= #

def make__sample():

    x_, y_, z_ = 0, 1, 2
    i_, s_, f_ = 3, 4, 5
    radius     = 0.8

    # ------------------------------------------------- #
    # --- [1] make coordinate                       --- #
    # ------------------------------------------------- #
    import nkUtilities.equiSpaceGrid as esg
    x1MinMaxNum = [ -1.0, 1.0, 101 ]
    x2MinMaxNum = [ -1.0, 1.0, 101 ]
    x3MinMaxNum = [  0.0, 0.0,   1 ]
    ret         = esg.equiSpaceGrid( x1MinMaxNum=x1MinMaxNum, x2MinMaxNum=x2MinMaxNum, \
                                     x3MinMaxNum=x3MinMaxNum, returnType = "structured" )

    Data                = np.zeros( (x3MinMaxNum[2],x2MinMaxNum[2],x1MinMaxNum[2],6) )
    Data[:,:,:,x_:z_+1] = ret
    radii               = np.sqrt( Data[:,:,:,x_]**2 + Data[:,:,:,y_]**2 )
    Data[:,:,:,     z_] = 1.0 - radii**2
    Data[:,:,:,     i_] = np.copy( Data[:,:,:,z_] )
    Data[:,:,:,     s_] = np.copy( Data[:,:,:,z_] )
    Data[:,:,:,     f_] = np.where( radii <= radius, 1, 0 )
    ( Data[...,z_] )[ np.where( radii > radius ) ] = 0.0
    
    # ------------------------------------------------- #
    # --- [2] save result                           --- #
    # ------------------------------------------------- #
    outFile   = "dat/input.dat"
    import nkUtilities.save__pointFile as spf
    spf.save__pointFile( outFile=outFile, Data=Data )

    # ------------------------------------------------- #
    # --- [3] show color map                        --- #
    # ------------------------------------------------- #
    import nkUtilities.cMapTri        as cmt
    pngFile = "png/input.png"
    Data_   = np.reshape( Data, (-1,6) )
    cmt.cMapTri( xAxis=Data_[:,x_], yAxis=Data_[:,y_], cMap=Data_[:,z_], \
                 pngFile=pngFile )

    # vtkFile = "png/input.vtp"
    # import nkVTKRoutines.convert__vtkPolySurface as vps
    # vps.convert__vtkPolySurface( Data=Data_, outFile=vtkFile )
    
    return()


# ========================================================= #
# ===   実行部                                          === #
# ========================================================= #

if ( __name__=="__main__" ):
    make__sample()
