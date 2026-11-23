import sys,os,math,struct,collections; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__))); import geo,old,ec
def plane_inv(p,P):
  o,U,V,_,_=geo.plane_params(p); d=[P[k]-o[k] for k in range(3)]
  dot=lambda a,b:sum(x*y for x,y in zip(a,b))
  aa,ab,bb=dot(U,U),dot(U,V),dot(V,V); r1,r2=dot(d,U),dot(d,V); det=aa*bb-ab*ab
  return ((bb*r1-ab*r2)/det,(aa*r2-ab*r1)/det)
c=collections.Counter(); shown=0
for l in open('cgmindex.tsv').read().splitlines():
  h=l.split('\t')[0]; s=open('cgm2/'+h+'.cgm','rb').read()
  if not old.old_records(s): continue
  o=geo.parse(s); C=geo.Ctx(o)
  for i,(t,refs,p,_) in o.items():
    if t!=14: continue
    e=refs[0]
    if o[e][0] not in (191,192,193): continue
    q=C.point14(i)[0][1]; m,_,_=ec.parse_ec(o[e][1],o[e][2]); P=None; tgt=None
    for cc,mp in zip(o[e][1],m):
      tc,cr,cp,_=o[cc]
      if not cr: continue
      if tc in (150,151):
        sid,f=C.pcurve(cc); S=C.surf(sid)
        if S: P=S(*f(mp(q)))
      if tc==160 and o[cr[0]][0]==200: tgt=(cc,cr[0],mp(q))
    if P is None or tgt is None: continue
    uv=plane_inv(o[tgt[1]][2],P); cp=o[tgt[0]][2]
    hits=[]
    for k in range(0,len(cp)-7):
      x=struct.unpack_from('<d',cp,k)[0]
      if abs(x-uv[0])<1e-6*max(1,abs(uv[0])): hits.append(('u',k))
      if abs(x-uv[1])<1e-6*max(1,abs(uv[1])): hits.append(('v',k))
    c[tuple(hits[:4])]+=1
    if shown<4: shown+=1; print(h,len(cp),'q',tgt[2],'uv',uv,hits)
print(c.most_common(10))
