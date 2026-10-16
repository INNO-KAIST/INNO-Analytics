#!/usr/bin/env python3
"""
d^n ligand-field (Tanabe-Sugano) calculator, built from first principles.

Hamiltonian:  H = sum_{i<j} e^2/r_ij  (Racah A,B,C)  +  V_oct(Dq)   [no spin-orbit]
Basis: all Slater determinants of d^n (complex spherical-harmonic orbitals m=-2..2, spin up/down).
Coulomb matrix elements: <ab|1/r|cd> = delta(ms_a,ms_c) delta(ms_b,ms_d) delta(m_a+m_b,m_c+m_d)
                         * sum_k c^k(l m_a, l m_c) c^k(l m_d, l m_b) F^k
   with Condon-Shortley c^k (Gaunt) coefficients; F^0=A+49F4.., Racah: A=F0-49F4, B=F2-5F4, C=35F4
   (F2 = F^2/49, F4 = F^4/441).
Octahedral crystal field (eg above t2g), complex basis:
   <+-2|V|+-2> = Dq, <+-1|V|+-1> = -4Dq, <0|V|0> = 6Dq, <2|V|-2> = <-2|V|2> = 5Dq
   => e_g = +6Dq, t2g = -4Dq, 10Dq = Delta_o.
States are classified by S (from <S^2>) and by the O-group irrep using the character of the C4(z)
rotation acting on orbital part only (chi(C4): A1=1, A2=-1, E=0, T1=1, T2=-1) and the orbital degeneracy.
Validated against textbook free-ion term energies (e.g. d3: 4P-4F = 15B, 2G-4F = 4B+3C) and the closed-form
d3 quartet energies 4T2 = 10Dq, 4T1(F,P) = 7.5B + 15Dq -/+ 0.5*sqrt(225B^2 - 180 B Dq + 100 Dq^2).
Tetrahedral: use Dq -> -(4/9)Dq_oct convention (i.e. pass negative Dq), labels then refer to Td (drop g).
"""
import numpy as np, itertools
from math import factorial, sqrt

def _w3j(j1,j2,j3,m1,m2,m3):
    if m1+m2+m3!=0 or abs(m1)>j1 or abs(m2)>j2 or abs(m3)>j3: return 0.0
    if j3<abs(j1-j2) or j3>j1+j2: return 0.0
    t1=j2-m1-j3; t2=j1+m2-j3; t3=j1+j2-j3; t4=j1-m1; t5=j2+m2
    tmin=max(0,t1,t2); tmax=min(t3,t4,t5)
    s=0.0
    for t in range(tmin,tmax+1):
        s+=(-1)**t/(factorial(t)*factorial(t-t1)*factorial(t-t2)*factorial(t3-t)*factorial(t4-t)*factorial(t5-t))
    pre=sqrt(factorial(j1+j2-j3)*factorial(j1-j2+j3)*factorial(-j1+j2+j3)/factorial(j1+j2+j3+1)*
             factorial(j1+m1)*factorial(j1-m1)*factorial(j2+m2)*factorial(j2-m2)*factorial(j3+m3)*factorial(j3-m3))
    return (-1)**(j1-j2-m3)*pre*s

def ck(k,l,m,lp,mp):
    # c^k(lm, l'm') = (-1)^m sqrt((2l+1)(2l'+1)) (l k l'; 0 0 0)(l k l'; -m m-m' m')
    return (-1)**m*sqrt((2*l+1)*(2*lp+1))*_w3j(l,k,lp,0,0,0)*_w3j(l,k,lp,-m,m-mp,mp)

ORB=[(m,s) for m in range(-2,3) for s in (+1,-1)]   # spin-orbitals (m, 2ms)

def _coul(a,b,c,d,F):
    (ma,sa),(mb,sb),(mc,sc),(md,sd)=ORB[a],ORB[b],ORB[c],ORB[d]
    if sa!=sc or sb!=sd or ma+mb!=mc+md: return 0.0
    return sum(ck(k,2,ma,2,mc)*ck(k,2,md,2,mb)*F[k] for k in (0,2,4))

def _cf(m,mp,Dq):
    if m==mp: return {2:Dq,-2:Dq,1:-4*Dq,-1:-4*Dq,0:6*Dq}[m]
    if abs(m-mp)==4: return 5*Dq
    return 0.0

class DN:
    def __init__(self,n):
        self.n=n
        self.dets=list(itertools.combinations(range(10),n))
        self.index={d:i for i,d in enumerate(self.dets)}
        N=len(self.dets)
        self.U=np.zeros((5,10))
        # precompute operator matrices per unit parameter
        self.Hk={k:np.zeros((N,N)) for k in (0,2,4)}
        self.Hcf=np.zeros((N,N))
        self.S2=np.zeros((N,N))
        for i,D in enumerate(self.dets):
            occ=list(D)
            # one-body CF
            for p in occ:
                for q in range(10):
                    if ORB[q][1]!=ORB[p][1]: continue
                    v=_cf(ORB[q][0],ORB[p][0],1.0)
                    if v==0: continue
                    new,ph=self._exc1(occ,p,q)
                    if new is not None: self.Hcf[self.index[new],i]+=ph*v
            # two-body Coulomb: sum_{p<q, r<s} <rs||pq> a_r+ a_s+ a_q a_p
            for p,q in itertools.combinations(occ,2):
                for r,s in itertools.combinations(range(10),2):
                    for k in (0,2,4):
                        F={0:0,2:0,4:0}; F[k]=1.0
                        v=_coul(r,s,p,q,F)-_coul(r,s,q,p,F)
                        if v==0: continue
                        new,ph=self._exc2(occ,p,q,r,s)
                        if new is not None: self.Hk[k][self.index[new],i]+=ph*v
            # S^2 = Sz^2 + Sz + S- S+   ... use S^2 = S-S+ + Sz(Sz+1)
            sz=sum(ORB[p][1] for p in occ)/2
            self.S2[i,i]+=sz*(sz+1)
            for p in occ:          # S+ : down -> up on same m, then S-
                if ORB[p][1]!=-1: continue
                q=self._partner(p)
                new,ph=self._exc1(occ,p,q)
                if new is None: continue
                occ2=list(new)
                for p2 in occ2:
                    if ORB[p2][1]!=+1: continue
                    q2=self._partner(p2)
                    new2,ph2=self._exc1(occ2,p2,q2)
                    if new2 is not None: self.S2[self.index[new2],i]+=ph*ph2
        # orbital angular momentum L^2 = (L+)^dag L+ + Lz^2 + Lz
        Lp=np.zeros((N,N)); Lz=np.zeros(N)
        for i,D in enumerate(self.dets):
            Lz[i]=sum(ORB[p][0] for p in D)
            for p in D:
                m,sp=ORB[p]
                if m==2: continue
                q=ORB.index((m+1,sp))
                new,ph=self._exc1(list(D),p,q)
                if new is not None: Lp[self.index[new],i]+=ph*sqrt(6-m*(m+1))
        self.L2=Lp.T@Lp+np.diag(Lz**2+Lz)
        self.C4=np.array([np.exp(-1j*np.pi/2*sum(ORB[p][0] for p in D)) for D in self.dets])
    @staticmethod
    def _partner(p): return p^1
    @staticmethod
    def _exc1(occ,p,q):
        if q!=p and q in occ: return None,0
        lst=list(occ); i=lst.index(p); ph=(-1)**i; lst.pop(i)
        j=sum(1 for x in lst if x<q); ph*=(-1)**j; lst.insert(j,q)
        return tuple(lst),ph
    @staticmethod
    def _exc2(occ,p,q,r,s):
        lst=list(occ); ph=1
        for x in (p,q):
            i=lst.index(x); ph*=(-1)**i; lst.pop(i)
        for x in (s,r):
            if x in lst: return None,0
            j=sum(1 for y in lst if y<x); ph*=(-1)**j; lst.insert(j,x)
        return tuple(lst),ph
    def H(self,B,C,Dq,A=0.0,alpha=0.0):
        F4=C/35.0; F2=B+5*F4; F0=A+49*F4
        # F^k (Slater) from normalized: F^2=49F2, F^4=441F4
        H=F0*self.Hk[0]+49*F2*self.Hk[2]+441*F4*self.Hk[4]+Dq*self.Hcf
        return H
    def free_ion_terms(self,B,C,alpha=0.0):
        """free-ion (Dq=0) terms: list of (E, 2S+1, L) sorted; accidental degeneracies lifted by tiny S^2/L^2 perturbation.
        alpha: Trees correction alpha*L(L+1) added (exact, since L^2 commutes with H at Dq=0)."""
        H=self.H(B,C,0.0)+alpha*self.L2
        w,v=np.linalg.eigh(H+1e-3*self.S2+1e-4*self.L2)
        out=[];i=0
        while i<len(w):
            j=i
            while j<len(w) and abs(w[j]-w[i])<1e-6: j+=1
            V=v[:,i:j]
            s2=np.trace(V.T@self.S2@V)/(j-i); l2=np.trace(V.T@self.L2@V)/(j-i)
            S=(-1+np.sqrt(1+4*s2))/2; L=(-1+np.sqrt(1+4*l2))/2
            e=np.trace(V.T@H@V)/(j-i)
            out.append((e,int(round(2*S+1)),int(round(L)),j-i)); i=j
        e0=min(o[0] for o in out)
        return sorted([(o[0]-e0,o[1],o[2]) for o in out])
    def _rotation_ops(self):
        if hasattr(self,'_R'): return self._R
        from scipy.linalg import expm
        ms=np.arange(-2,3)
        Lz=np.diag(ms).astype(complex)
        Lp=np.zeros((5,5),complex)
        for i,m in enumerate(ms[:-1]): Lp[i+1,i]=sqrt(6-m*(m+1))
        Lx=(Lp+Lp.T)/2; Ly=(Lp-Lp.T)/(2j)
        def D(axis,theta):
            n=np.array(axis,float); n/=np.linalg.norm(n)
            return expm(-1j*theta*(n[0]*Lx+n[1]*Ly+n[2]*Lz))
        ops={'C3':D((1,1,1),2*np.pi/3),'C2':D((0,0,1),np.pi),'C4':D((0,0,1),np.pi/2),"C2'":D((1,1,0),np.pi)}
        self._R={}
        for k,d in ops.items():
            U=np.zeros((10,10),complex)
            for a,(ma,sa) in enumerate(ORB):
                for b,(mb,sb) in enumerate(ORB):
                    if sa==sb: U[a,b]=d[ma+2,mb+2]
            N=len(self.dets); R=np.zeros((N,N),complex)
            for j,Dj in enumerate(self.dets):
                for i,Di in enumerate(self.dets):
                    R[i,j]=np.linalg.det(U[np.ix_(Di,Dj)])
            self._R[k]=R
        return self._R
    def levels(self,B,C,Dq,tol=1e-6):
        """return list of (E rel. to ground, 2S+1, O-irrep, total degeneracy), accidental degeneracies decomposed by characters"""
        R=self._rotation_ops()
        CT={'A1':[1,1,1,1,1],'A2':[1,1,1,-1,-1],'E':[2,-1,2,0,0],'T1':[3,0,-1,1,-1],'T2':[3,0,-1,-1,1]}
        cls=[1,8,3,6,6]
        w,v=np.linalg.eigh(self.H(B,C,Dq))
        out=[];i=0
        while i<len(w):
            j=i
            while j<len(w) and abs(w[j]-w[i])<1e-5*max(1,abs(w[i])): j+=1
            V=v[:,i:j]
            s2=V.T@self.S2@V; sv,su=np.linalg.eigh(s2)
            k=0
            while k<len(sv):
                l=k
                while l<len(sv) and abs(sv[l]-sv[k])<1e-6: l+=1
                W=V@su[:,k:l]; S=(-1+np.sqrt(1+4*sv[k]))/2; mult=int(round(2*S+1))
                chi=[ (l-k)/mult ]+[np.real(np.trace(W.conj().T@R[c]@W))/mult for c in ('C3','C2','C4',"C2'")]
                for ir,ct in CT.items():
                    nrep=int(round(sum(g*a*b for g,a,b in zip(cls,ct,chi))/24))
                    for _ in range(nrep): out.append([w[i],mult,ir,mult*ct[0]])
                k=l
            i=j
        E0=out[0][0]
        return [(round(e-E0,3),m,ir,d) for e,m,ir,d in out]
    def _levels_old(self,B,C,Dq,tol=1e-6):
        """return list of (E, 2S+1, irrep, degeneracy) sorted, E relative to lowest"""
        w,v=np.linalg.eigh(self.H(B,C,Dq))
        out=[]; i=0
        while i<len(w):
            j=i
            while j<len(w) and abs(w[j]-w[i])<1e-6*max(1,abs(w[i]))+1e-7: j+=1
            V=v[:,i:j]
            s2=np.trace(V.T@self.S2@V)/(j-i); S=(-1+np.sqrt(1+4*s2))/2; mult=int(round(2*S+1))
            g=(j-i)//mult if mult else 0
            chi=np.real(np.trace(V.conj().T@np.diag(self.C4)@V))/mult
            irrep={1:{1:'A1',-1:'A2'},2:{0:'E'},3:{1:'T1',-1:'T2'}}.get(g,{}).get(int(round(chi)),f'?g{g}chi{chi:.2f}')
            out.append([w[i],mult,irrep,j-i])
            i=j
        # accidental degeneracies (e.g. Dq=0 free ion or exact crossings) can merge irreps; caller should avoid exact crossings
        E0=out[0][0]
        return [(round(e-E0,3),m,ir,d) for e,m,ir,d in out]

def d3_quartets(B,Dq):
    r=0.5*np.sqrt(225*B**2-180*B*Dq+100*Dq**2)
    return {'4T2':10*Dq,'4T1(F)':7.5*B+15*Dq-r,'4T1(P)':7.5*B+15*Dq+r}

if __name__=='__main__':
    d3=DN(3)
    # free-ion check (Dq tiny to lift nothing: use Dq=0 and group)
    B,C=918.,3850.
    w=np.linalg.eigvalsh(d3.H(B,C,0.0)); w=np.round(w-w[0],3)
    print('d3 free-ion distinct term energies/B:',sorted(set(np.round(w/B,4))))
    print('expected 4P=15, 2G=4+3C/B=',4+3*C/B)
    L=d3.levels(B,C,2000.)
    print(L[:8]); print(d3_quartets(B,2000.))
