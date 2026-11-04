"""T=0 window: is the tetrahedron the ground state at zH4 = c*zG ? (J=1)"""
import numpy as np
from numpy.polynomial import legendre as L
from scipy.optimize import minimize
rng=np.random.default_rng(0)
def kap(K,zG,zH4,J=1.0):
    return np.array([0,J-0.6*zG,2/3*K-4/7*zH4,-0.4*zG,-8/35*zH4])
def Ebody(x,c):
    S=x.reshape(4,3);S=S/np.linalg.norm(S,axis=1)[:,None]
    X=S@S.T; iu=np.triu_indices(4,1)
    return L.legval(X[iu],c).sum()
def is_T(K,zG,zH4,nstart=40):
    c=kap(K,zG,zH4); eT=6*L.legval(-1/3,c); best=np.inf
    for _ in range(nstart):
        r=minimize(Ebody,rng.normal(size=12),args=(c,),method='BFGS')
        best=min(best,r.fun)
    return best>eT-1e-6, eT, best
if __name__=="__main__":
    import sys
    cfac=float(sys.argv[1]) if len(sys.argv)>1 else 9.70
    Ks=[0.5,1,2,3,4,6]; zGs=[0.02,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.8]
    print("c =",cfac," rows K, cols zG:",zGs)
    for K in Ks:
        row=[]
        for zG in zGs:
            ok,_,_=is_T(K,zG,cfac*zG); row.append('T' if ok else '.')
        # analytic bounds
        tc=0.9*K/(cfac-0.3); tf=(1+2*K/3)/(7/9+20/27*cfac)
        print(f"K={K:>3}: {' '.join(row)}   analytic zG < min({tc:.3f}, {tf:.3f})")
