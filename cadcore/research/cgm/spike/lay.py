import sys,os,collections,struct; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__))); import old,vi
def refs(b):
  c,k,_=vi.rd(b,0); out=[]
  for _ in range(c): x,k,_=vi.rd(b,k); out.append(x)
  return out,k
FAM=lambda v: 'A' if ('V5R1' in v and any(x in v for x in ['V5R17-JMB','V5R17-JMC','V5R18','V5R19'])) or 'R20' in v else 'B'
if __name__=='__main__':
  T=int(sys.argv[1]); c=collections.Counter()
  for f,s,r,e in old.files():
    v=old.version(s); fam=FAM(v)
    for (o,fl,t,l,i,b) in r:
      if t==T: rr,k=refs(b); c[(fam,len(rr),l-k)]+=1
  for k,n in sorted(c.items()): print(n,k)
