import numpy as np, window as w, twobody as t
def zgmax(K,c,lo=0.005,hi=1.0):
    if not w.is_T(K,lo,c*lo)[0]: return 0
    for _ in range(14):
        m=0.5*(lo+hi)
        if w.is_T(K,m,c*m,nstart=30)[0]: lo=m
        else: hi=m
    return lo
def reff(K,zG,zH4):
    G=zG/4;H=zH4/4;base=t.frame(np.ones(3));ref=t.body(base,0,1);d=[]
    for sig in (1,-1):
        for a in np.linspace(0,2*np.pi/3,12,endpoint=False):
            Sb=t.body(base,a,sig)
            d.append((sig,a,t.energy(ref,Sb,ref,1,K,G,H,Sb),t.relax_corr(ref,Sb,ref,1,K,G,H,Sb)))
    d=np.array(d);X=t.fit(d,(0,3,6))
    ce=np.linalg.lstsq(X,d[:,2]+d[:,3],rcond=None)[0]
    J4=-ce[2];J2=-ce[3];J1=-ce[1]
    return J4/J2, J1/J2, J4/(512/729*H)
for K in [0.5,1,2,3,4,6,8]:
    zm=zgmax(K,9.70)
    zG=0.5*zm
    r,j12,_=reff(K,zG,9.70*zG)
    print(f"K={K}: zG_max={zm:.3f} (T_x={1.768*zm:.3f} J at edge) | at zG={zG:.3f}: r_eff={r:.3f} (rigid 4.311), J1/J2={j12:.3f} (rigid 1.25)")
