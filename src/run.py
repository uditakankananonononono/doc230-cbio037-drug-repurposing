import gzip,json,numpy as np,pandas as pd,torch,torch.nn as nn,scipy.sparse as sp
from sklearn.decomposition import TruncatedSVD
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score,average_precision_score
torch.set_num_threads(2)
E=pd.read_csv('/tmp/het/edges.sif.gz',sep='\t'); N=pd.read_csv('/tmp/het/nodes.tsv',sep='\t')
comp=sorted(N[N.kind=='Compound'].id); dis=sorted(N[N.kind=='Disease'].id); gene=sorted(N[N.kind=='Gene'].id); se=sorted(N[N.kind=='Side Effect'].id)
ci={x:i for i,x in enumerate(comp)}; di={x:i for i,x in enumerate(dis)}; gi={x:i for i,x in enumerate(gene)}; si={x:i for i,x in enumerate(se)}
def mat(me,a,b,A,B,swap=False):
    d=E[E.metaedge==me]; r=[A[x] for x in d.source]; c=[B[x] for x in d.target]
    return sp.csr_matrix((np.ones(len(r)),(r,c)),shape=(len(A),len(B)))
CG=sum(mat(m,0,0,ci,gi) for m in ('CbG','CuG','CdG')); CG.data[:]=1
d=E[E.metaedge.isin(['DaG','DuG','DdG'])]
DG=sp.csr_matrix((np.ones(len(d)),([di[x] for x in d.source],[gi[x] for x in d.target])),shape=(len(dis),len(gene))); DG.data[:]=1
CS=mat('CcSE',0,0,ci,si); CS.data[:]=1
def pairs(me):
    d=E[E.metaedge==me]; return set((ci[a],di[b]) for a,b in zip(d.source,d.target))
CT=pairs('CtD'); CP=pairs('CpD'); ALL=CT|CP
print('CtD',len(CT),'CpD',len(CP),flush=True)
svd=TruncatedSVD(64,random_state=0).fit(sp.vstack([CG,DG])); CZ=svd.transform(CG); DZ=svd.transform(DG)
CSZ=TruncatedSVD(32,random_state=0).fit_transform(CS)
sc=lambda Z:(Z-Z.mean(0))/(Z.std(0)+1e-9)
CZ,DZ,CSZ=sc(CZ),sc(DZ),sc(CSZ)
def feats(P,mode='both'):
    a=np.array([p[0] for p in P]); b=np.array([p[1] for p in P])
    parts=[]
    if mode in('both','gene'): parts+= [CZ[a],DZ[b],CZ[a]*DZ[b]]
    if mode in('both','se'): parts+= [CSZ[a]]
    if mode=='se': parts+=[DZ[b]]
    return np.hstack(parts).astype(np.float32)
cos=lambda a,b:(CZ[a]*DZ[b]).sum(1)/(np.linalg.norm(CZ[a],axis=1)*np.linalg.norm(DZ[b],axis=1)+1e-9)
def sample_neg(rng,cs,ds,n,excl,w=None):
    out=set(); cs=np.array(cs); ds=np.array(ds)
    pc=None if w is None else w[0]; pd_=None if w is None else w[1]
    while len(out)<n:
        c=rng.choice(cs,size=n,p=pc); dd=rng.choice(ds,size=n,p=pd_)
        for x in zip(c,dd):
            if x not in excl and x not in out: out.add(x)
            if len(out)>=n: break
    return list(out)
def mlp(X,y,seed,epochs=30):
    torch.manual_seed(seed); net=nn.Sequential(nn.Linear(X.shape[1],128),nn.ReLU(),nn.Dropout(.3),nn.Linear(128,128),nn.ReLU(),nn.Dropout(.3),nn.Linear(128,1))
    opt=torch.optim.Adam(net.parameters(),1e-3,weight_decay=1e-4); Xt=torch.tensor(X); yt=torch.tensor(y,dtype=torch.float32); lf=nn.BCEWithLogitsLoss()
    for ep in range(epochs):
        net.train(); pm=torch.randperm(len(Xt))
        for i in range(0,len(Xt),64):
            b=pm[i:i+64]; opt.zero_grad(); lf(net(Xt[b]).squeeze(-1),yt[b]).backward(); opt.step()
    net.eval(); return lambda X:net(torch.tensor(X)).squeeze(-1).detach().numpy()
def split(kind,rng,target):
    pos=sorted(target); cdeg=np.zeros(len(comp)); ddeg=np.zeros(len(dis))
    cwith=sorted(set(p[0] for p in pos)); dwith=sorted(set(p[1] for p in pos))
    if kind=='random':
        idx=rng.permutation(len(pos)); k=int(.8*len(pos)); tr=[pos[i] for i in idx[:k]]; te=[pos[i] for i in idx[k:]]
        trc=list(range(len(comp))); trd=list(range(len(dis))); tec=trc; ted=trd
    elif kind=='cold_compound':
        hc=set(rng.choice(cwith,size=int(.2*len(cwith)),replace=False)); tr=[p for p in pos if p[0] not in hc]; te=[p for p in pos if p[0] in hc]
        trc=[c for c in range(len(comp)) if c not in hc]; trd=list(range(len(dis))); tec=sorted(hc); ted=trd
    else:
        hd=set(rng.choice(dwith,size=int(.2*len(dwith)),replace=False)); tr=[p for p in pos if p[1] not in hd]; te=[p for p in pos if p[1] in hd]
        trc=list(range(len(comp))); trd=[x for x in range(len(dis)) if x not in hd]; tec=trc; ted=sorted(hd)
    return tr,te,trc,trd,tec,ted
def run(kind,seed,target=CT,hard=False,shuffle=False,mode='both',res_extra=False):
    rng=np.random.default_rng(seed); tr,te,trc,trd,tec,ted=split(kind,rng,target)
    w_tr=w_te=None
    if hard:
        cd=np.bincount([p[0] for p in tr],minlength=len(comp))+1.; dd=np.bincount([p[1] for p in tr],minlength=len(dis))+1.
        w_tr=(cd[trc]/cd[trc].sum(),dd[trd]/dd[trd].sum()); w_te=(cd[tec]/cd[tec].sum(),dd[ted]/dd[ted].sum())
    ntr=sample_neg(rng,trc,trd,len(tr),ALL,w_tr); nte=sample_neg(rng,tec,ted,len(te),ALL,w_te)
    P=tr+ntr; y=np.array([1]*len(tr)+[0]*len(ntr)); 
    if shuffle: y=rng.permutation(y)
    Q=te+nte; yq=np.array([1]*len(te)+[0]*len(nte))
    X=feats(P,mode); Xq=feats(Q,mode); s={}
    cd=np.bincount([p[0] for p in tr],minlength=len(comp)); dd=np.bincount([p[1] for p in tr],minlength=len(dis))
    s['B1']=np.array([(cd[a]+1)*(dd[b]+1) for a,b in Q],float)
    s['B3']=cos(np.array([q[0] for q in Q]),np.array([q[1] for q in Q]))
    lr=LogisticRegression(max_iter=300).fit(X,y); s['B2']=lr.decision_function(Xq)
    f=mlp(X,y,seed); s['M']=f(Xq)
    out={k:dict(auroc=roc_auc_score(yq,v),auprc=average_precision_score(yq,v)) for k,v in s.items()}
    out['n_test_pos']=len(te); out['n_train_pos']=len(tr)
    if res_extra: out['_tr']=(tr,te,trc,tec,ted,f,cd)
    return out
def agg(rs):
    ks=[k for k in rs[0] if k in('B1','B2','B3','M')]
    return {k:{m:[round(float(np.mean([r[k][m] for r in rs])),3),round(float(np.std([r[k][m] for r in rs])),3)] for m in('auroc','auprc')} for k in ks}
R={}; SEEDS=range(5)
R['d1']=dict(n_compounds=len(comp),n_diseases=len(dis),n_genes=len(gene),CtD=len(CT),CpD=len(CP),compounds_with_CtD=len(set(p[0] for p in CT)),diseases_with_CtD=len(set(p[1] for p in CT)))
for name,kw in [('d2_d9_cold_compound',dict(kind='cold_compound')),('d3_random',dict(kind='random')),('d4_cold_disease',dict(kind='cold_disease')),
  ('d5_shuffle_control',dict(kind='cold_compound',shuffle=True)),('d6_degree_matched_neg_cold_compound',dict(kind='cold_compound',hard=True)),
  ('d6_degree_matched_neg_random',dict(kind='random',hard=True)),('d7_gene_only',dict(kind='cold_compound',mode='gene')),('d7_se_only',dict(kind='cold_compound',mode='se')),
  ('d8_CpD_cold_compound',dict(kind='cold_compound',target=CP))]:
    rs=[run(seed=s,**kw) for s in SEEDS]; R[name]=agg(rs); print(name,R[name]['M'],R[name]['B1'],R[name]['B3'],flush=True)
bb=max(R['d2_d9_cold_compound'][k]['auroc'][0] for k in('B1','B2','B3'))
R['d2_gate_M_beats_best_baseline_by_0.02']=bool(R['d2_d9_cold_compound']['M']['auroc'][0]>=bb+0.02)
R['d3_gate_gap_le_0.05']=bool(R['d3_random']['M']['auroc'][0]-R['d2_d9_cold_compound']['M']['auroc'][0]<=0.05)
# D10 hits@10 per held-out disease, cold compound
hits={'M':[],'B1':[],'B3':[]}
for s in SEEDS:
    o=run('cold_compound',s,res_extra=True); tr,te,trc,tec,ted,f,cd=o['_tr']
    rng=np.random.default_rng(100+s); nocT=[c for c in range(len(comp)) if c not in set(p[0] for p in CT)]
    cand=sorted(set(tec)|set(rng.choice(nocT,200,replace=False)))
    for d_ in sorted(set(p[1] for p in te)):
        tp=set(p[0] for p in te if p[1]==d_); P=[(c,d_) for c in cand]; X=feats(P)
        sc_={'M':f(X),'B3':cos(np.array(cand),np.full(len(cand),d_)),'B1':cd[np.array(cand)].astype(float)+1}
        for k,v in sc_.items():
            top=np.array(cand)[np.argsort(-v+1e-9*np.random.default_rng(0).random(len(v)))[:10]]; hits[k].append(int(len(tp&set(top))>0))
R['d10_hits_at_10_frac_diseases']={k:[round(float(np.mean(v)),3),len(v)] for k,v in hits.items()}
R['d10_note']='candidate set = held-out CtD compounds + 200 random compounds with no CtD; rankings are benchmark scores, not validated repurposing claims'
json.dump(R,open('results/results.json','w'),indent=1); print(json.dumps(R,indent=1))
