import sys,os,collections,struct; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__))); import old,lay
T=int(sys.argv[1]); N=int(sys.argv[2]) if len(sys.argv)>2 else 6; want=sys.argv[3] if len(sys.argv)>3 else None
seen=collections.Counter()
for f,s,r,e in old.files():
  v=old.version(s).strip()[:9]
  if want and want not in v: continue
  for (o,fl,t,l,i,b) in r:
    if t==T:
      rr,k=lay.refs(b); key=l-k
      if seen[key]>=N: continue
      seen[key]+=1; p=b[k:]
      nd=len(p)//8
      print(v,i,rr,key,'|',' '.join('%.6g'%x for x in struct.unpack_from('<%dd'%nd,p,0)) if (len(p)%8==0 and 'hex' not in sys.argv) else p.hex(' '))
