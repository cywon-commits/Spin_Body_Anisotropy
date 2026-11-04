"""Two-body checks for the C3/S4 effective models (rigid and relaxed).
Square lattice z=4. Body a at origin, b at +x; their other 3 neighbours each are
rigid tetrahedra with the same axis at alpha=0, sigma=+1.
Intra: J x + K x^2 ; inter (per bond, all 16 pairs): -(G x^3 + H4 x^4), G = zG/4.
"""
import numpy as np, sys
V = np.array([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]])/np.sqrt(3)

def frame(axis):
    z = axis/np.linalg.norm(axis)
    x = np.cross(z,[0.3,0.7,0.1]); x/=np.linalg.norm(x); y=np.cross(z,x)
    return V@np.array([x,y,z]).T          # body with 'axis' along lab z

def Rz(a):
    c,s=np.cos(a),np.sin(a); return np.array([[c,-s,0],[s,c,0],[0,0,1]])

def body(base, alpha, sigma):
    return sigma*(base@Rz(alpha).T)

def intra(S,J,K):
    e=0.
    for i in range(4):
        for j in range(i+1,4):
            x=S[i]@S[j]; e+=J*x+K*x*x
    return e
def inter(A,B,G,H):
    X=A@B.T; return -np.sum(G*X**3+H*X**4)

def energy(Sa,Sb,ref,J,K,G,H,refb=None):
    refb=ref if refb is None else refb
    return (intra(Sa,J,K)+intra(Sb,J,K)+inter(Sa,Sb,G,H)
            +3*inter(Sa,ref,G,H)+3*inter(Sb,refb,G,H))

def tangent(S):
    T=[]
    for s in S:
        e1=np.cross(s,[0.31,0.52,0.79]); e1/=np.linalg.norm(e1); e2=np.cross(s,e1)
        T.append((e1,e2))
    return T

def relax_corr(Sa,Sb,ref,J,K,G,H,refb,h=1e-4):
    S0=np.vstack([Sa,Sb]); T=tangent(S0); n=16
    def E(u):
        S=S0.copy()
        for k in range(8):
            v=S0[k]+u[2*k]*T[k][0]+u[2*k+1]*T[k][1]; S[k]=v/np.linalg.norm(v)
        return energy(S[:4],S[4:],ref,J,K,G,H,refb)
    g=np.zeros(n); Hm=np.zeros((n,n)); e0=E(np.zeros(n))
    for i in range(n):
        ei=np.zeros(n); ei[i]=h
        g[i]=(E(ei)-E(-ei))/(2*h)
        for j in range(i,n):
            ej=np.zeros(n); ej[j]=h
            Hm[i,j]=Hm[j,i]=(E(ei+ej)-E(ei-ej)-E(-ei+ej)+E(-ei-ej))/(4*h*h)
    # rigid-rotation modes of a and b
    Q=[]
    for b0 in (0,4):
        for ax in np.eye(3):
            v=np.zeros(n)
            for k in range(4):
                d=np.cross(ax,S0[b0+k]); kk=b0+k
                v[2*kk]=d@T[kk][0]; v[2*kk+1]=d@T[kk][1]
            Q.append(v)
    Q=np.linalg.qr(np.array(Q).T)[0]; P=np.eye(n)-Q@Q.T
    Hp=P@Hm@P; gp=P@g
    w,U=np.linalg.eigh(Hp); keep=w>1e-6
    if np.any(w< -1e-6): print("WARNING negative mode",w[w<-1e-6])
    return -0.5*sum((U[:,k]@gp)**2/w[k] for k in range(n) if keep[k])

def run(axis,zG,zH4,J=1.,K=0.5,relax=True):
    G=zG/4; H=zH4/4; base=frame(axis); ref=body(base,0,1)
    out=[]
    for sig in (1,-1):
        for a in np.linspace(0,2*np.pi/3 if axis is not None else 1,25,endpoint=False):
            Sa=body(base,0,1); Sb=body(base,a,sig)
            e=energy(Sa,Sb,ref,J,K,G,H)-energy(Sa,ref,ref,J,K,G,H)  # bond excess
            c=relax_corr(Sa,Sb,ref,J,K,G,H) if relax else 0.
            out.append((sig,a,e,c))
    return np.array(out)

def fit(d, harmonics):
    """fit y(sig,a) = sum_k [A_k + B_k sig] cos(k a)"""
    rows=[];y=[]
    for s,a,e,c in d:
        rows.append([f(s)*np.cos(k*a) for k in harmonics for f in (lambda s:1,lambda s:s)])
    return np.array(rows)

if __name__=="__main__":
    for zG,zH4 in [(1.0,0.0),(1.0,1.0)]:
        G=zG/4;H=zH4/4
        for name,axis,hs in [("C3",np.ones(3),(0,3,6)),("S4",np.array([0,0,1.]),(0,2,4,6,8))]:
            base=frame(axis); ref=body(base,0,1)
            d=[]
            per=2*np.pi/3 if name=="C3" else np.pi
            for sig in (1,-1):
                for a in np.linspace(0,per,24,endpoint=False):
                    Sa=body(base,0,1); Sb=body(base,a,sig)
                    e=energy(Sa,Sb,ref,1,0.5,G,H,Sb)
                    c=relax_corr(Sa,Sb,ref,1,0.5,G,H,Sb)
                    d.append((sig,a,e,c))
            d=np.array(d); X=fit(d,hs)
            ce=np.linalg.lstsq(X,d[:,2],rcond=None)[0]; cc=np.linalg.lstsq(X,d[:,3],rcond=None)[0]
            rese=np.abs(X@ce-d[:,2]).max(); resc=np.abs(X@cc-d[:,3]).max()
            print(f"\nzG={zG} zH4={zH4} {name}: rigid fit resid {rese:.1e}, relax fit resid {resc:.1e}")
            for i,k in enumerate(hs):
                print(f"  cos({k}a): rigid const {ce[2*i]:+.5f} sig {ce[2*i+1]:+.5f} | relax const {cc[2*i]:+.6f} sig {cc[2*i+1]:+.6f}")
            if name=="C3":
                print(f"  predicted rigid per bond: J1={160/81*G:.5f} J2={128/81*G:.5f} J4={512/729*H:.5f}")
            else:
                print(f"  predicted rigid: cos2 {-32/9*G:.5f}, cos4 {-32/81*H:.5f}")
