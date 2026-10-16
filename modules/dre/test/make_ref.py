import numpy as np, json
from scipy.integrate import solve_ivp
from pa3 import params_lee, rhs, all_steady
S={'No2':(2.25,0.35,0.6),'No3':(8,3.25,6.3),'No4':(8,0.25,0.5),'No5':(8,0.0,0.0),'No6':(20,6.0,14.0),'No7':(20,3.6,7.2)}
I_cm=np.logspace(3,5.1,40)   # W/cm2
out={}
for k,(c,f2,f3) in S.items():
    p=params_lee(c,f2,f3)
    st=[]
    for I in I_cm:
        # 독립 방법: 바닥상태에서 Radau 로 충분히 오래 적분 (정상상태 도달 확인)
        sol=solve_ivp(lambda t,n: rhs(n,p,I*1e4),[0,300.0],[1,0,0],method='LSODA',rtol=1e-10,atol=1e-20)
        a=sol.y[:,-1]; b=solve_ivp(lambda t,n: rhs(n,p,I*1e4),[0,600.0],a,method='LSODA',rtol=1e-10,atol=1e-20).y[:,-1]
        assert abs(b[2]/a[2]-1)<1e-7, (k,I)
        roots=all_steady(p,I*1e4,ngrid=6000)
        st.append([[b[1],b[2]]]+[[x,y] for x,y in roots])
    # 동역학: 3개 파워
    kin={}
    for I in [3e3,7e3,2e4]:
        sol=solve_ivp(lambda t,n: rhs(n,p,I*1e4),[0,2.0],[1,0,0],method='Radau',rtol=1e-10,atol=1e-18,dense_output=True)
        t=np.logspace(-5,np.log10(2.0),60); kin[str(I)]={'t':t.tolist(),'n2':sol.sol(t)[1].tolist(),'n3':sol.sol(t)[2].tolist()}
    out[k]={'c':c,'w2nr':p['w2nr'],'w3nr':p['w3nr'],'I':I_cm.tolist(),'steady':st,'kin':kin}
json.dump(out,open('lee_ref.json','w'))
print({k:max(len(s) for s in v['steady']) for k,v in out.items()})
