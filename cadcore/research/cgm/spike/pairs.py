import json,os,sys,collections,pickle; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__))); import cpt,geo,old
CC={5:13,6:14,24:150,25:151,32:160,35:191,36:192,37:193,39:200,40:201,42:203,43:204,46:212,52:231,55:241,91:500,92:502,93:510,94:511,95:512,96:513,97:520,98:521,100:524,101:525}
idx=[l.split('\t') for l in open('cgmindex.tsv').read().splitlines()]
byname=collections.defaultdict(list)
for h,_,p in idx: byname[os.path.basename(p).lower()].append(h)
P=collections.defaultdict(list); stats=collections.Counter()
for n,hs in byname.items():
  olds=[];cps=[]
  for h in hs:
    s=open('cgm2/'+h+'.cgm','rb').read()
    if old.old_records(s):
      try: olds.append((h,geo.parse(s)))
      except Exception: pass
    else:
      x=cpt.records(s)
      if x: cps.append((h,{i:(t,fl,b) for (i,t,fl,b) in x[0]},x[1][2]))
  for hc,C,L in cps:
    for ho,O in olds:
      common=set(C)&set(O)
      agree=sum(1 for i in common if CC.get(C[i][0])==O[i][0])
      if common and agree>0.9*len(common) and len(common)>0.8*max(len(C),len(O)):
        stats[('paired',L)]+=1
        for i in common:
          if CC.get(C[i][0])==O[i][0]: P[C[i][0]].append((n,i,O[i],C[i][2],L))
        break
pickle.dump(dict(P),open('pairs.pkl','wb'))
print(stats, {t:len(v) for t,v in sorted(P.items())})
