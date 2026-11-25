"""pairs old-encoding records with compact records of the same part saved by another release, by
aligning the class sequences (compact ids are renumbered from 2012 on)"""
import json,os,sys,collections,pickle,difflib; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__))); import cpt,geo,old
from pairs import CC
idx=[l.split('\t') for l in open('cgmindex.tsv').read().splitlines()]
byname=collections.defaultdict(list)
for h,_,p in idx: byname[os.path.basename(p).lower()].append(h)
P=collections.defaultdict(list); stats=collections.Counter()
for n,hs in byname.items():
  olds=[];cps=[]
  for h in hs:
    s=open('cgm2/'+h+'.cgm','rb').read()
    if old.old_records(s):
      try: o=geo.parse(s); olds.append((h,[(i,)+o[i] for i in sorted(o)]))
      except Exception: pass
    else:
      x=cpt.records(s)
      if x and (x[1][2]<79 or x[1][2]==81): cps.append((h,x[0],x[1][2]))
  for hc,R,L in cps:
    best=None
    for ho,O in olds:
      a=[t for (i,t,r,p,_) in O]; b=[CC.get(t,-t) for (i,t,fl,bb) in R]
      sm=difflib.SequenceMatcher(None,a,b,autojunk=False)
      m=sum(bl.size for bl in sm.get_matching_blocks())
      if m>0.9*max(len(a),len(b)) and (best is None or m>best[0]): best=(m,O,sm)
    if not best: continue
    stats[('paired',L)]+=1
    m,O,sm=best
    for bl in sm.get_matching_blocks():
      for j in range(bl.size):
        oi=O[bl.a+j]; ci=R[bl.b+j]
        P[ci[1]].append((n,oi[0],(oi[1],oi[2],oi[3],oi[4]),ci[3],L))
pickle.dump(dict(P),open('pairs2.pkl','wb'))
print(stats, {t:len(v) for t,v in sorted(P.items())})
