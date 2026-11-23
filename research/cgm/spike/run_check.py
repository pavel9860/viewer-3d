import sys,os,collections; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__))); import geo,old
bo=[l.split('\t')[0] for l in open('brep_only.tsv').read().splitlines()]
tot=collections.Counter(); W=[]
for h in bo:
  try: s=open('cgm2/'+h+'.cgm','rb').read()
  except: continue
  if not old.old_records(s): continue
  st,w=geo.check(geo.parse(s)); tot.update(st); W+=[(h,)+x for x in w[:1]]
for k,v in sorted(tot.items(),key=lambda x:-x[1])[:int(sys.argv[1]) if len(sys.argv)>1 else 30]: print(v,k)
for x in W[:6]: print(x)
