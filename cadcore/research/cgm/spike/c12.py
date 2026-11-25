"""2012+ (08 b5 75) compact -> canonical records, with cross-record derivations (surface boxes,
p-curve ranges, edge-curve ranges). Returns objs {id: (old class, refs, payload, 0)} for geo.check2."""
import struct,cpt,c11,plchk,collections
HYP={"b2":"o0","b3":"sym"}
D=c11.D; h=c11.h; I=c11.I; V=c11.V
class Ctx:
  def __init__(s,R):
    s.R={i:(t,b) for (i,t,fl,b) in R}; s.objs={}; s.err=collections.Counter(); s.box={}; s.prange={}; s.erange={}; s.ptleft=collections.Counter()
  def run(s):
    order=lambda t: 0 if t in (39,40,41,42,43,44,45,46,48,52,55,56,59,2,4,14) else 1 if t in (24,25,32) else 2 if t in (35,36,37) else 3
    for i,(t,b) in sorted(s.R.items(),key=lambda x:order(x[1][0])):
      try: s.objs[i]=s.conv(i,t,b)
      except KeyError: s.err[('nokey',t)]+=1
      except Exception as e: s.err[(t,str(e)[:40])]+=1
    return s.objs
  def conv(s,i,t,b):
    if t==39: return s.plane(i,b)
    if t==24: return s.line(i,b)
    if t in (35,36,37): return s.ec(i,t,b)
    if t==6: return s.point(i,b)
    if t==97: return s.shell(i,b)
    ot,refs,p=c11.conv(t,b,i,True,False)
    if ot==201:
      v=struct.unpack_from('<17d',p); s.box[i]=v[10:14]
    if ot==151:
      v=struct.unpack_from('<2d2i5d',p); s.prange[i]=(v[5],v[6])
    if ot==160:
      s.prange[i]=struct.unpack_from('<2d',p,len(p)-16)
    return (ot,refs,p,0)
  def plane(s,i,b):
    o,U,V_,sc,box=plchk.plane(b)
    s.box[i]=box
    return (200,[],D(*o,*U,*V_,*sc,*box),0)
  def line(s,i,b):
    r=c11.R(b,i); sid=r.refs(1)[0]; f=r.u8()
    box=s.box.get(sid)
    def need():
      if box is None: raise ValueError('line: no box for surface %d'%s.R[sid][0])
      return box
    H=HYP
    if f&3:   # +v / +u: [c, t0, t1]; bit 4: c = box start (cross axis); bit 2: t0 = box start;
              # bit 3: t1 = box end (t0 derived) or -t0
      vert=f&1; d=(0.,1.) if vert else (1.,0.)
      if f&0x20: raise ValueError('line bit5 positive %02x'%f)
      c=(need()[0] if vert else need()[2]) if f&16 else r.f()
      t0=(need()[2] if vert else need()[0]) if f&4 else r.f()
      if f&8: t1=(need()[3] if vert else need()[1]) if f&4 else -t0
      else: t1=r.f()
      o=(c,0.) if vert else (0.,c)
    elif f&0xc0:
      vert=f&64; d=(0.,-1.) if vert else (-1.,0.)
      if f&0x30:
        bx=need(); lo,hi=(bx[2],bx[3]) if vert else (bx[0],bx[1]); t0,t1=-hi,-lo
      else:
        t0=r.f()
        if f&8:
          if H['b3']=='sym': t1=-t0
          else: t0,t1=0.0,t0
        else: t1=r.f()
      o=(0.,0.) if f&4 else tuple(r.fs(2))
    else:
      if f&0x30: raise ValueError('line general bits %02x'%f)
      t0=r.f()
      if f&8:
        if H['b3']=='sym': t1=-t0
        else: t0,t1=0.0,t0
      else: t1=r.f()
      o=(0.,0.) if f&4 else tuple(r.fs(2)); d=tuple(r.fs(2))
    if r.left(): raise ValueError('line %d left'%r.left())
    s.prange[i]=(t0,t1)
    return (150,[sid],D(*o,*d,t0,t1),0)
  def ec(s,i,t,b):
    r=c11.R(b,i); refs=r.refs(); hd=r.u8(); x=r.u8()
    if x>>2!=0x21: raise ValueError('ec x %02x'%x)
    n=1 if hd&8 else r.count()
    p=V(n); tail=False; first=None
    for j in range(n):
      m=0 if hd&4 else r.sint()
      last=j==n-1
      if j==0 and m in (0,1) and n==1:
        pr=s.prange.get(refs[0])
        if x&1 or x&2:
          if pr is None: raise ValueError('ec: no p-curve range (%d)'%s.R[refs[0]][0])
        a0=pr[0] if x&1 else r.f(); a1=pr[1] if x&2 else r.f()
      else:
        a0,a1=r.f(),r.f()
      if first is None: first=(a0,a1)
      if m==0:
        tol=0.0 if (last and hd&1) else r.f(); p+=h(m)+D(a0,a1,tol); continue
      tol=r.f(); p+=h(m)+D(a0,a1,tol)
      if m==1:
        b0,b1=r.f(),r.f(); tb=0.0 if (last and hd&1 and not tail) else r.f(); p+=D(b0,b1,tb); continue
      tail=True; p+=D(r.f(),r.f())
      if m==3: p+=D(r.f(),r.f())
      elif m==5: N=r.count(); p+=V(N)+D(*r.fs(2*N))
      elif m!=2: raise ValueError('ec map %d'%m)
    if tail:
      c0=first[0] if x&1 else r.f(); c1=first[1] if x&2 else r.f(); ct=0.0 if hd&1 else r.f(); p+=D(c0,c1,ct)
    if r.left(): raise ValueError('ec %d left'%r.left())
    s.erange[i]=first
    return (c11.CLS[t],refs,p+b'\0\0\0\0',0)
  def point(s,i,b):
    # refs; code byte: bit 7: refs 0 and 1 share the parameter; bits 2-3: refs[0] (2 start, 1 end of
    # the edge-curve range, 0 stored f64 + int unless the group's curve is a line or circle); bits 0-1
    # and 4-5: the other refs (not decoded: their parameters follow from the edge curve's maps)
    r=c11.R(b,i); refs=r.refs(); cd=r.u8()
    code=cd>>2&3
    grp=refs[:2] if cd&0x80 else refs[:1]
    if code==0:
      q=r.f()
      k=1 if s.R[grp[-1]][0] in (24,25,36) else r.sint()
    else:
      er=s.erange.get(refs[0]) or s.prange.get(refs[0])
      if er is None: raise ValueError('point: no range for %d'%s.R[refs[0]][0])
      q=er[0] if code==2 else er[1]; k=1
    if r.left()==0 or True: pass
    p=V(1)+D(q)+I(k)
    s.ptleft[(cd>>4&3, cd&3, len(refs), r.left()==0)]+=1
    return (14,[refs[0]],p,0)
  def shell(s,i,b):
    try: ot,refs,p=c11.conv(97,b,i,True,False); return (ot,refs,p,0)
    except Exception:
      r=c11.R(b,i); refs=r.refs()
      if b[r.k]!=0: raise
      b2=b[:r.k]+b[r.k+1:]
      ot,refs,p=c11.conv(97,b2,i,True,False); return (ot,refs,p,0)
