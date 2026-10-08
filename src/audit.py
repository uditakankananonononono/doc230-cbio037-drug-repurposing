import json
SRC=open('src/run.py').read(); exec(SRC.split("R={}; SEEDS")[0])
G4=np.asarray((CG>0).sum(1)).ravel().astype(float)
def run2(seed,matched):
    rng=np.random.default_rng(seed); tr,te,trc,trd,tec,ted=split('cold_compound',rng,CT)
    if matched:
        ctr=sorted(set(p[0] for p in tr)); cte=sorted(set(p[0] for p in te)); ntr=sample_neg(rng,ctr,trd,len(tr),ALL); nte=sample_neg(rng,cte,ted,len(te),ALL)
    else:
        ntr=sample_neg(rng,trc,trd,len(tr),ALL); nte=sample_neg(rng,tec,ted,len(te),ALL)
    P=tr+ntr; y=np.array([1]*len(tr)+[0]*len(ntr)); Q=te+nte; yq=np.array([1]*len(te)+[0]*len(nte))
    X=feats(P); Xq=feats(Q); s={}
    s['B4']=np.array([G4[a] for a,b in Q]); s['B3']=cos(np.array([q[0] for q in Q]),np.array([q[1] for q in Q]))
    s['B2']=LogisticRegression(max_iter=300).fit(X,y).decision_function(Xq); s['M']=mlp(X,y,seed)(Xq)
    return {k:float(roc_auc_score(yq,v)) for k,v in s.items()}
R={}
for nm,mt in (('original_random_negatives',False),('compound_matched_negatives',True)):
    rs=[run2(s,mt) for s in range(5)]; R[nm]={k:dict(mean=float(np.mean([r[k] for r in rs])),sd=float(np.std([r[k] for r in rs]))) for k in rs[0]}; print(nm,json.dumps(R[nm]),flush=True)
R['Q1_B4_ge_M_minus_0.02']=bool(R['original_random_negatives']['B4']['mean']>=R['original_random_negatives']['M']['mean']-0.02)
R['Q2_M_gt_0.60']=bool(R['compound_matched_negatives']['M']['mean']>0.60)
hits={k:[] for k in ('M','B3','B4','random')}
for s in range(5):
    rng=np.random.default_rng(s); tr,te,trc,trd,tec,ted=split('cold_compound',rng,CT)
    # refit M exactly as original run (neg sampling stream identical to run(): split, then negatives)
    ntr=sample_neg(rng,trc,trd,len(tr),ALL,None); nte=sample_neg(rng,tec,ted,len(te),ALL,None)
    P=tr+ntr; y=np.array([1]*len(tr)+[0]*len(ntr)); f=mlp(feats(P),y,s)
    rng2=np.random.default_rng(100+s); nocT=[c for c in range(len(comp)) if c not in set(p[0] for p in CT)]
    cand=np.array(sorted(set(tec)|set(rng2.choice(nocT,200,replace=False))))
    for d_ in sorted(set(p[1] for p in te)):
        tp=set(p[0] for p in te if p[1]==d_); X=feats([(c,d_) for c in cand])
        for k,v in (('M',f(X)),('B3',cos(cand,np.full(len(cand),d_))),('B4',G4[cand])):
            top=cand[np.argsort(-v+1e-9*np.random.default_rng(0).random(len(v)))[:10]]; hits[k].append(int(len(tp&set(top))>0))
        hits['random'].append(float(np.mean([len(tp&set(np.random.default_rng(i).choice(cand,10,replace=False)))>0 for i in range(20)])))
R['Q3_hits_at_10']={k:[float(np.mean(v)),len(v)] for k,v in hits.items()}
R['Q3_B4_ge_M_minus_0.05']=bool(R['Q3_hits_at_10']['B4'][0]>=R['Q3_hits_at_10']['M'][0]-0.05)
json.dump(R,open('results/audit.json','w'),indent=1); print(json.dumps(R,indent=1)); print('done')
