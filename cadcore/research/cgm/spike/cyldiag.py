import sys,os,math,collections; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__))); import geo,old,ec
bo=[l.split('\t')[0] for l in open('brep_only.tsv').read().splitlines()]
c=collections.Counter(); ex={}
for h in bo:
  try: s=open('cgm2/'+h+'.cgm','rb').read()
  except: continue
  if not old.old_records(s): continue
  o=geo.parse(s); C=geo.Ctx(o)
  for i,(t,refs,p,_) in o.items():
    if t!=14: continue
    e=refs[0]
    if o[e][0] not in (191,192,193): continue
    q=C.point14(i)[0][1]; m,_,_=ec.parse_ec(o[e][1],o[e][2])
    P=None; cyl=None
    for cc,mp in zip(o[e][1],m):
      tc,cr,cp,_=o[cc]
      if not cr: continue
      st=o[cr[0]][0]
      if tc in (150,151) and st==200: sid,f=C.pcurve(cc); P=C.surf(sid)(*f(mp(q)))
      if tc==150 and st==201: cyl=(cr[0],C.pcurve(cc)[1](mp(q)))
    if P is None or cyl is None: continue
    v=geo.D(o[cyl[0]][2],0,len(o[cyl[0]][2])//8); O,X,Y,r=v[0:3],v[3:6],v[6:9],v[9]
    d=[P[k]-O[k] for k in range(3)]; th=math.atan2(sum(d[k]*Y[k] for k in range(3)),sum(d[k]*X[k] for k in range(3)))
    u=cyl[1][0]; off=v[16] if len(v)>16 else 0
    diff=(th-(u-off)/(r*v[14]))%(2*math.pi)
    key=round(diff/(math.pi/2),3)%4
    k2=(key, round((v[10]-off)/(r*math.pi/2),3), len(v), v[14])
    c[k2]+=1; ex.setdefault(k2,(h,i))
for k,n in sorted(c.items(),key=lambda x:-x[1])[:25]: print(n,k,ex[k])
