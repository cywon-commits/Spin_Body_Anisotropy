import numpy as np
from numba import njit
from sbaniso.c3model import wolff_sigma_phi, sweep, measure, seed, work_arrays, new_state, NOBS
@njit
def run(s,p,L,beta,c,stack,inc,mem,n,mode,ncl,out):
    o=np.zeros(10)
    for it in range(n):
        if mode==0:
            tot=0
            while tot<L*L:
                tot+=wolff_sigma_phi(s,p,L,beta,c,stack,inc,mem)
        elif mode==1:
            for k in range(ncl):
                wolff_sigma_phi(s,p,L,beta,c,stack,inc,mem)
        else:
            sweep(s,p,L,beta,c,1.0)
            tot=0
            while tot<L*L:
                tot+=wolff_sigma_phi(s,p,L,beta,c,stack,inc,mem)
        measure(s,p,L,c,o); out[it,0]=o[0]/L**2; out[it,1]=(o[1]/L**2)**2
L=6; c=np.array([1.0,0,0]); b=1/2.5
for mode,name in [(0,"until-N"),(1,"fixed 3 clusters"),(2,"metro+until-N")]:
    s,p=new_state(L,np.random.default_rng(1)); st,inc,mem=work_arrays(L); seed(2)
    out=np.zeros((200000,2)); run(s,p,L,b,c,st,inc,mem,200000,mode,3,out)
    bl=out.reshape(50,-1,2).mean(1); print(name, bl.mean(0).round(4), (bl.std(0)/7).round(4))
print("exact  [-1.2801  0.5286]")
