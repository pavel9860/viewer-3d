import sys,os,struct,math,collections
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__))); import old,vi
def parse(s):
  r,e=old.old_records(s)
  objs={}
  for (o,fl,t,l,i,b) in r:
    c,k,_=vi.rd(b,0); refs=[]
    for _ in range(c): x,k,_=vi.rd(b,k); refs.append(x)
    objs[i]=(t,refs,b[k:],fl)
  return objs
D=lambda p,k,n=1: struct.unpack_from('<%dd'%n,p,k)
def add(a,b): return [x+y for x,y in zip(a,b)]
def mul(a,s): return [x*s for x in a]
def plane_params(p):
  if len(p)==120: v=D(p,0,15); return v[0:3],v[3:6],v[6:9],v[9:11],v[11:15]
  if len(p)%8==1:
    o=D(p,0,3); m=p[24]; n=(len(p)-25)//8; vals=list(D(p,25,n)); comps=[]
    for j in range(6): comps.append(0.0 if m>>j&1 else vals.pop(0))
    return o,comps[0:3],comps[3:6],vals[0:2],vals[2:6]
  raise ValueError('plane len %d'%len(p))
class Ctx:
  def __init__(s,objs): s.o=objs
  def surf(s,i):
    t,refs,p,_=s.o[i]
    if t==200:
      o,U,V,sc,box=plane_params(p)
      return lambda u,v: add(o,add(mul(U,u),mul(V,v)))
    if t==201:
      v=D(p,0,len(p)//8); o,X,Y=v[0:3],v[3:6],v[6:9]; r=v[9]
      Z=[X[1]*Y[2]-X[2]*Y[1],X[2]*Y[0]-X[0]*Y[2],X[0]*Y[1]-X[1]*Y[0]]
      def f(u,w):
        th=u/(r*v[14])
        return add(o,add(mul(X,r*math.cos(th)),add(mul(Y,r*math.sin(th)),mul(Z,w))))
      return f
    if t==241:
      return s.surf(refs[0])
    if t==203:
      v=D(p,0,len(p)//8); o,A,B,Cz=v[0:3],v[3:6],v[6:9],v[9:12]; rs=v[12]
      def f(u,w):
        u,w=u/rs,w/rs
        cu,su,cw,sw=math.cos(u),math.sin(u),math.cos(w),math.sin(w)
        return [o[k]+cw*(cu*A[k]+su*B[k])+sw*Cz[k] for k in range(3)]
      return f
    if t==204:
      v=D(p,0,len(p)//8); o,X,Y,Z=v[0:3],v[3:6],v[6:9],v[9:12]; R,r=v[12],v[13]
      def f(u,w):
        u,w=u/R,w/r
        cu,su,cw,sw=math.cos(u),math.sin(u),math.cos(w),math.sin(w); q=R+r*cw
        return [o[k]+q*(cu*X[k]+su*Y[k])+r*sw*Z[k] for k in range(3)]
      return f
    if t==212:
      return s.surf(refs[0])
    if t==231:
      import nurbs
      S,_=nurbs.surface(p)
      return lambda u,v: nurbs.surf_eval(S,u,v)
    return None
  def pcurve(s,i):
    t,refs,p,_=s.o[i]
    if t==150:
      if len(p)==48: v=D(p,0,6); o,d=v[0:2],v[2:4]
      else: v=D(p,0,8); o,d=v[0:2],v[3:5]
      return refs[0],lambda q:(o[0]+q*d[0],o[1]+q*d[1])
    if t==151:
      c=D(p,0,2); v=D(p,24,(len(p)-24)//8); r=v[0]
      sc,off=v[3],v[4] if len(v)>4 else 0
      return refs[0],lambda q:(c[0]+r*math.cos(q/r/sc),c[1]+r*math.sin(q/r/sc))
    if t==160:
      import herm
      H,_=herm.parse_herm(p)
      return refs[0],lambda q:tuple(herm.herm_eval(H,q))
    return None
  def point14(s,i):
    t,refs,p,_=s.o[i]
    k=0
    if (len(p)-1)%12==0: n,k,_=vi.rd(p,0)
    else: n=struct.unpack_from('<I',p,0)[0]; k=4
    out=[]
    for j in range(n): out.append((refs[j],D(p,k+12*j)[0],struct.unpack_from('<i',p,k+12*j+8)[0]))
    return out
def check(objs,verbose=False):
  return check2(objs)
def check2(objs):
  C=Ctx(objs); st=collections.Counter(); worst=[]
  for i,(t,refs,p,_) in objs.items():
    if t!=14: continue
    ec=refs[0]
    if objs[ec][0] not in (192,193,191): st['first %d'%objs[ec][0]]+=1; continue
    q=C.point14(i)[0][1]
    pts=[]
    import ec as ECm
    try: maps,_,rest=ECm.parse_ec(objs[ec][1],objs[ec][2])
    except Exception as ex: st['ecfail']+=1; continue
    if False: st['ecrest %d'%len(rest)]+=1
    for c,mp in zip(objs[ec][1],maps):
      qq=mp(q)
      tc=objs[c][0]
      pc=C.pcurve(c) if tc in (150,151,160) else None
      if not pc: st['skip %d'%tc]+=1; continue
      sid,f=pc; S=C.surf(sid)
      if not S: st['nosurf %d'%objs[sid][0]]+=1; continue
      pts.append((tc,objs[sid][0],S(*f(qq))))
    if len(pts)>=2:
      d=max(math.dist(a[2],b[2]) for a in pts for b in pts)
      key=(objs[ec][0],)+tuple(sorted(set((a[0],a[1]) for a in pts)))
      st[('ok' if d<1e-3 else 'bad',key)]+=1
      if d>=1e-3 and len(worst)<5: worst.append((i,d,pts))
  return st,worst
def check_old(objs):
  C=Ctx(objs); st=collections.Counter(); worst=[]
  for i,(t,refs,p,_) in objs.items():
    if t!=14: continue
    pts=[]
    for (c,q,flag) in C.point14(i):
      tc=objs[c][0]
      pc=C.pcurve(c) if tc in (150,151,160) else None
      if not pc: st['skip %d'%tc]+=1; continue
      sid,f=pc; S=C.surf(sid)
      if not S: st['nosurf %d'%objs[sid][0]]+=1; continue
      uv=f(q); pts.append((tc,objs[sid][0],S(*uv)))
    if len(pts)>=2:
      d=max(math.dist(a[2],b[2]) for a in pts for b in pts)
      key=tuple(sorted(set((a[0],a[1]) for a in pts)))
      st[('ok' if d<1e-3 else 'bad',key)]+=1
      if d>=1e-3 and len(worst)<5: worst.append((i,d,pts))
  return st,worst
if __name__=='__main__':
  s=open(sys.argv[1],'rb').read(); st,w=check(parse(s))
  for k,v in sorted(st.items(),key=str): print(v,k)
  for x in w: print(x)
