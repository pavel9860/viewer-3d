import sys,os,math,collections; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__))); import geo,old,ec
T=int(sys.argv[1]); N=int(sys.argv[2]) if len(sys.argv)>2 else 12
out=[]
for l in open('cgmindex.tsv').read().splitlines():
  h=l.split('\t')[0]; s=open('cgm2/'+h+'.cgm','rb').read()
  if not old.old_records(s): continue
  o=geo.parse(s); C=geo.Ctx(o)
  def carrier(sid):
    while o[sid][0] in (241,212): sid=o[sid][1][0]
    return sid
  for i,(t,refs,p,_) in o.items():
    if t!=14: continue
    e=refs[0]
    if o[e][0] not in (191,192,193): continue
    q=C.point14(i)[0][1]; m,_,_=ec.parse_ec(o[e][1],o[e][2]); P=None; tgt=None
    for cc,mp in zip(o[e][1],m):
      tc,cr,cp,_=o[cc]
      if not cr or tc not in (150,151,160): continue
      sid=cr[0]; car=carrier(sid)
      if o[car][0]==T: tgt=(car,C.pcurve(cc)[1](mp(q)))
      elif o[car][0] in (200,201,231):
        S=C.surf(sid)
        if S: P=S(*C.pcurve(cc)[1](mp(q)))
    if P and tgt: out.append((h,tgt[0],tgt[1],P,geo.D(o[tgt[0]][2],0,len(o[tgt[0]][2])//8)))
    if len(out)>=N: break
  if len(out)>=N: break
for h,sid,uv,P,v in out:
  o_,A,B,Cz=v[0:3],v[3:6],v[6:9],v[9:12]
  d=[P[k]-o_[k] for k in range(3)]; dot=lambda a,b:sum(x*y for x,y in zip(a,b))
  x,y,z=dot(d,A)/dot(A,A),dot(d,B)/dot(B,B),dot(d,Cz)/dot(Cz,Cz)
  print(h,sid,'uv',[round(t,4) for t in uv],'local',[round(t,4) for t in (x,y,z)],'atan2',round(math.atan2(y,x),4),'params',[round(t,4) for t in v[12:]])
