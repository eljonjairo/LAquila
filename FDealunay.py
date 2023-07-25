import numpy as np
from scipy.spatial import Delaunay

def Facetx(X,Y,Z,xb):
    
    ix = np.where(X == xb)
    ye = Y[ix]
    ze = Z[ix]
    YZ  = np.array([ye,ze]).transpose()
    tri = Delaunay(YZ).simplices
    ntri = int(tri.size/3)
    triX = np.zeros((ntri,3))
    jx = np.asarray(ix).transpose()
   
    for itri in range(0,ntri):
        xtri = tri[itri,0]
        ytri = tri[itri,1]
        ztri = tri[itri,2]
        triX[itri,0] = jx[xtri]
        triX[itri,1] = jx[ytri]
        triX[itri,2] = jx[ztri]
         
    return triX

def Facety(X,Y,Z,yb):
    
    iy = np.where(Y == yb)
    xe = X[iy]
    ze = Z[iy]
    XZ  = np.array([xe,ze]).transpose()
    tri = Delaunay(XZ).simplices
    ntri = int(tri.size/3)
    triY = np.zeros((ntri,3))
    jy = np.asarray(iy).transpose()
    
    for itri in range(0,ntri):
        xtri = tri[itri,0]
        ytri = tri[itri,1]
        ztri = tri[itri,2]
        triY[itri,0] = jy[xtri]
        triY[itri,1] = jy[ytri]
        triY[itri,2] = jy[ztri]
      
    return triY

def Facetz(X,Y,Z,zb):
    
    iz = np.where(Z == zb)
    xe = X[iz]
    ye = Y[iz]
    XY  = np.array([xe,ye]).transpose()
    tri = Delaunay(XY).simplices
    ntri = int(tri.size/3)
    triZ = np.zeros((ntri,3))
    jz = np.asarray(iz).transpose()
    
    for itri in range(0,ntri):
        xtri = tri[itri,0]
        ytri = tri[itri,1]
        ztri = tri[itri,2]
        triZ[itri,0] = jz[xtri]
        triZ[itri,1] = jz[ytri]
        triZ[itri,2] = jz[ztri]
      
    return triZ

