"""content pairing: a compact record is the old record of the same class whose 8-byte doubles it
contains (most votes, unique); topology records then pair through matched references"""
import os,sys,collections,pickle,struct; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__))); import cpt,geo,old
from pairs import CC
NAMES=['471510001_hdmi','672982090_usb','sj1-3555ng','c-5558065-1-a1-3d','pin_strip-2.54','473460001_microusb','5009980900_sd_reader']
idx=[l.split('\t') for l in open('cgmindex.tsv').read().splitlines()]
P=collections.defaultdict(list); st=collections.Counter()
for nm in NAMES:
  hs=[h for h,_,p in idx if os.path.basename(p).lower()==nm+'.catpart']
  O=None; C=None
  for h in hs:
    s=open('cgm2/'+h+'.cgm','rb').read()
    if old.old_records(s): O=geo.parse(s)
    else:
      x=cpt.records(s)
      if x and x[1][2] in (74,75,76,77,78,81): C=x[0]
  if O is None or C is None: st[('missing',nm)]+=1; continue
  # index: (class, double bytes) -> old ids
  ix=collections.defaultdict(set)
  for i,(t,r,p,_) in O.items():
    for k in range(0,len(p)-7):
      d=p[k:k+8]
      if d.count(0)<6: ix[(t,d)].add(i)
  match={}
  for (ci,t,fl,b) in C:
    ot=CC.get(t)
    if ot is None: continue
    votes=collections.Counter()
    for k in range(0,len(b)-7):
      d=b[k:k+8]
      if d.count(0)<6:
        for oi in ix.get((ot,d),()): votes[oi]+=1
    if votes:
      (oi,v),*rest=votes.most_common(2)+[(None,0)]
      if v>rest[0][1]: match[ci]=oi
  st[(nm,'matched',len(match),'of',len(C))]+=1
  for (ci,t,fl,b) in C:
    if ci in match: P[t].append((nm,match[ci],O[match[ci]],b,ci))
pickle.dump(dict(P),open('pairs3.pkl','wb'))
for k,v in st.items(): print(k)
print({t:len(v) for t,v in sorted(P.items())})
