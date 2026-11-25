import collections,cpt,struct,new3
def vec3(z,unit,it):
  w=[0.0 if z>>j&1 else None for j in range(3)]
  if unit:
    free=[j for j in range(3) if w[j] is None]
    if len(free)!=1: raise ValueError('unit with %d free'%len(free))
    w[free[0]]=1.0
  return [x if x is not None else next(it) for x in w]
def plane(b):
  v,k=cpt.sint(b,0)
  n=(len(b)-k)//8
  if (len(b)-k)%8: raise ValueError('len')
  it=iter(struct.unpack_from('<%dd'%n,b,k))
  o=vec3(v&7,0,it); U=vec3(v>>4&7,v>>3&1,it); V=vec3(v>>8&7,v>>7&1,it)
  sc=[1.0,1.0] if v>>11&1 else [next(it),next(it)]
  box=[-1e7,1e7,-1e7,1e7] if v>>12&1 else [next(it) for _ in range(4)]
  rest=list(it)
  if rest: raise ValueError('rest %d'%len(rest))
  if v>>13: raise ValueError('high bits')
  return o,U,V,sc,box
if __name__=='__main__':
  c=collections.Counter(); ex={}
  for f,s,R in cpt.files():
    r,d,k=new3.chain(s); L=d[2]
    if L not in (74,75,76,77,78,81): continue
    for (i,t,fl,b) in R:
      if t!=39: continue
      try:
        o,U,V,sc,box=plane(b)
        on=abs(sum(x*x for x in U)-1)<1e-9 and abs(sum(x*x for x in V)-1)<1e-9 and abs(sum(a*b for a,b in zip(U,V)))<1e-9
        key=('ok' if on else 'notortho', tuple(sc)==(1,1), box==[-1e7,1e7,-1e7,1e7])
      except Exception as e: key=('err',str(e)[:30], b[0]&3)
      c[key]+=1; ex.setdefault(key,(f,i,b.hex()[:120]))
  for k,v in sorted(c.items(),key=str): print(v,k,ex[k])
