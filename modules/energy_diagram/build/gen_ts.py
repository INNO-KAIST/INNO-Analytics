import numpy as np, json, time
from dn_ligand_field import DN
DQ=[round(0.05*i,2) for i in range(1,81)]      # Dq/B 0.05..4.0
CB=[round(3.4+0.2*i,1) for i in range(12)]     # C/B 3.4..5.6
out={}
for n in (3,5,8):
    t=time.time(); d=DN(n); B=1000.0
    grid={}   # label -> [ci][qi]
    for ci,cb in enumerate(CB):
        for qi,dq in enumerate(DQ):
            H=d.H(B,cb*B,dq*B); E0=np.linalg.eigvalsh(H)[0]
            L=d.levels(B,cb*B,dq*B)
            cnt={}
            for e,m,ir,deg in sorted(L,key=lambda x:x[0]):
                k=f"{m}{ir}"; i=cnt.get(k,0); cnt[k]=i+1
                lab=f"{k}#{i}"
                grid.setdefault(lab,[[None]*len(DQ) for _ in CB])
                grid[lab][ci][qi]=round((e)/B,3)   # rel to ground
    keep={}
    for lab,g in grid.items():
        vals=[v for row in g for v in row if v is not None]
        if len(vals)<len(DQ)*len(CB): continue
        if min(vals)<=60: keep[lab]=g
    out[f"d{n}"]={"labels":sorted(keep,key=lambda l:keep[l][6][30]),"E":[keep[l] for l in sorted(keep,key=lambda l:keep[l][6][30])]}
    print(n,len(grid),len(keep),time.time()-t)
json.dump({"DQ":DQ,"CB":CB,"grids":out,"note":"E/B relative to ground at each (Dq/B, C/B); octahedral, no spin-orbit. Generated from dn_ligand_field.py (exact d^n diagonalization)."},open("ts_grids.json","w"),separators=(",",":"))
