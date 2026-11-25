"""2011Spring (levels 79/80) compact record -> canonical (old varint-dialect) record:
(old class, absolute refs, old payload bytes)."""
import struct,cpt,vi
CLS={2:10,4:12,5:13,6:14,14:100,24:150,25:151,26:152,32:160,35:191,36:192,37:193,39:200,40:201,41:202,42:203,43:204,
     44:210,45:211,46:212,48:214,52:231,55:241,91:500,92:502,93:510,94:511,95:512,96:513,97:520,98:521,100:524,101:525}
class R:
  def __init__(s,b,i,ci=True): s.b=b; s.k=0; s.i=i; s.ci=ci
  def i16(s):
    if s.ci: return s.sint()
    x=struct.unpack_from('<h',s.b,s.k)[0]; s.k+=2; return x
  def i32(s):
    if s.ci: return s.sint()
    x=struct.unpack_from('<i',s.b,s.k)[0]; s.k+=4; return x
  def sint(s): v,s.k=cpt.sint(s.b,s.k); return v
  def ref(s): return s.i-s.sint()
  def count(s): v,s.k,_=vi.rd(s.b,s.k); return v
  def refs(s,n=None):
    if n is None: n=s.count()
    return [s.ref() for _ in range(n)]
  def u8(s): x=s.b[s.k]; s.k+=1; return x
  def raw(s,n):
    x=s.b[s.k:s.k+n]; s.k+=n
    if len(x)!=n: raise ValueError('short')
    return x
  def f(s): x=struct.unpack_from('<d',s.b,s.k)[0]; s.k+=8; return x
  def fs(s,n): return [s.f() for _ in range(n)]
  def rest(s): x=s.b[s.k:]; s.k=len(s.b); return x
  def left(s): return len(s.b)-s.k
  def bits(s,n):
    nb=(n+7)//8; x=int.from_bytes(s.b[s.k:s.k+nb],'little'); s.k+=nb
    return [(x>>j)&1 for j in range(n)]
def V(n): # old varint count
  if n<128: return bytes([0x80|n])
  bs=n.to_bytes(4,'little'); m=0; out=b''
  for j in range(4):
    if bs[j]: m|=1<<j; out+=bs[j:j+1]
  return bytes([m<<3])+out
D=lambda *xs: struct.pack('<%dd'%len(xs),*xs)
h=lambda *xs: struct.pack('<%dh'%len(xs),*xs)
I=lambda *xs: struct.pack('<%di'%len(xs),*xs)
def pm(x): return 1 if x else -1
def conv(t,b,i,ci=True,dr=True):
  r=R(b,i,ci); ot=CLS.get(t)
  if ot is None: raise KeyError('code %d'%t)
  if ot==13: refs=r.refs(); p=b''
  elif ot==14:
    refs=r.refs(); p=V(len(refs))
    for _ in refs: q=r.f(); p+=D(q)+I(r.i32())
  elif ot==510: refs=r.refs(1); p=bytes([r.u8()])
  elif ot==511: refs=([0]+r.refs(4)) if dr else r.refs(5); p=bytes([r.u8()])
  elif ot==512: refs=r.refs(); p=h(r.i16())
  elif ot==520:
    refs=r.refs(); n=len(refs); bt=r.bits(4+n)
    kind={1:1,3:3,8:-1}.get(bt[0]|bt[1]<<1|bt[2]<<2|bt[3]<<3)
    if kind is None: raise ValueError('shell kind')
    p=V(n)+h(kind,*[pm(x) for x in bt[4:]])+h(r.i16())
  elif ot==521:
    refs=r.refs(); n=r.count()
    if 2*n+1!=len(refs): raise ValueError('loop refs %d n %d'%(len(refs),n))
    hd=r.u8(); ab={0x41:(1,1),0x50:(1,-1),0xc1:(-1,1),0xd0:(-1,-1)}.get(hd)
    if ab is None: raise ValueError('loop head %02x'%hd)
    bt=r.bits(3*n); vals=[-1]+[pm(x) for x in bt]   # x0 implied -1
    c=r.i16()
    p=V(n)+h(*ab)+h(*vals[:3*n])+h(vals[3*n],c)
    if c:   # block: 4 f64, 4 ints, u8, 24 bytes (not decoded)
      p+=r.raw(32)+I(*[r.i32() for _ in range(4)])+bytes([r.u8()])+r.raw(24)
  elif ot in (513,):refs=r.refs(); p=b''
  elif ot==525: refs=r.refs(); p=I(r.i32())+h(r.i16())
  elif ot==524: refs=r.refs(1); p=h(r.i16())
  elif ot==500: refs=r.refs(); p=r.rest()   # payload not used by the decoder
  elif ot==502: refs=r.refs(); p=r.rest()
  elif ot in (191,192,193):
    refs=r.refs(); hd=r.u8()
    n=1 if hd&8 else r.count()
    k0=r.k; sols=[]
    for opt in range(1<<n):     # bit j: map-1 entry j omits its identity range (b0 b1 tol)
      try:
        r.k=k0; p=V(n); tail=False
        for j in range(n):
          m=0 if hd&4 else r.i16()
          a0,a1=r.f(),r.f()
          last = j==n-1
          if m!=1 and opt>>j&1: raise ValueError('opt')
          if m==0:
            tol=0.0 if (last and hd&1) else r.f(); p+=h(m)+D(a0,a1,tol); continue
          tol=r.f(); p+=h(m)+D(a0,a1,tol)
          if m==1:
            if opt>>j&1: p+=D(a0,a1,tol); continue
            b0,b1=r.f(),r.f(); tb=0.0 if (last and hd&1 and not tail) else r.f()
            p+=D(b0,b1,tb); continue
          if opt>>j&1: raise ValueError('opt')
          tail=True
          b0,b1=r.f(),r.f(); p+=D(b0,b1)
          if m==3: p+=D(r.f(),r.f())
          elif m==5:
            N=r.count(); p+=V(N)+D(*r.fs(2*N))
          elif m!=2: raise ValueError('ec map %d'%m)
        if tail:
          c0,c1=r.f(),r.f(); ct=0.0 if hd&1 else r.f(); p+=D(c0,c1,ct)
        if r.left()==8: r.raw(8)   # two i32 (raw-int streams; not decoded)
        if r.left()==0: sols.append((r.k,p))
      except (ValueError,IndexError,struct.error): pass
    if len(sols)!=1: raise ValueError('ec layouts %d'%len(sols))
    r.k,p=sols[0]; p+=b'\0\0\0\0'
  elif ot==150:
    refs=r.refs(1); form=r.i16()
    if form==0: o=r.fs(2); d=r.fs(2); t0,t1=r.f(),r.f()
    elif form==1: o=[r.f(),0.0]; d=[0.0,1.0]; t0,t1=r.f(),r.f()
    elif form==2: o=[0.0,r.f()]; d=[1.0,0.0]; t0,t1=r.f(),r.f()
    else: raise ValueError('line form %d'%form)
    p=D(*o,*d,t0,t1)
  elif ot==151:
    refs=r.refs(1); c=r.fs(2)
    def fits(k):   # r t0 t1 at k, flag byte, [s], off
      return k+24<len(b) and b[k+24] in (0,1) and len(b)-k==25+(8 if b[k+24]==0 else 0)+8
    plain,ints=fits(r.k),fits(r.k+8)
    if plain and ints and b[r.k:r.k+8]==I(1,1): plain=False
    if plain==ints: raise ValueError('circle variant %d %d'%(plain,ints))
    if ints: r.raw(8)   # the variant with two i32 (1, 1)
    rad,t0,t1=r.fs(3); fl=r.u8()
    s=1.0 if fl==1 else r.f()
    if fl not in (0,1): raise ValueError('circle flag %d'%fl)
    off=r.f()
    p=D(*c)+I(1,1)+D(rad,t0,t1,s,off)
  elif ot in (10,100,202,203,204): refs=[]; p=r.rest()
  elif ot in (12,211): refs=r.refs(1); p=r.rest()
  elif ot==200:
    refs=[]; o=r.fs(3); m=r.u8()
    UV=[0.0 if m>>j&1 else r.f() for j in range(6)]
    sc=[1.0,1.0] if m&64 else r.fs(2)
    box=[-1e7,1e7,-1e7,1e7] if m&128 else r.fs(4)
    p=D(*o,*UV,*sc,*box)
  elif ot==201:
    refs=[]; o=r.fs(3); m=r.u8()
    XY=[0.0 if m>>j&1 else r.f() for j in range(6)]
    if m>>6: raise ValueError('cyl mask %02x'%m)
    v=r.fs(5); f2=r.u8()
    su=1.0 if f2&1 else r.f(); sv=1.0 if f2&2 else r.f(); x=v[1] if f2&4 else r.f()
    if f2>>3: raise ValueError('cyl tail %02x'%f2)
    p=D(*o,*XY,*v,su,sv,x)
  elif ot==160:
    refs=r.refs(1); deg=r.i32(); n=r.i32(); z=r.u8()
    if z!=4: raise ValueError('herm z %02x'%z)
    kn=r.fs(n); mm=r.u8()
    if mm==1: mu=[deg+1]+[(deg+1)//2]*(n-2)+[deg+1]
    else: raise ValueError('herm mults %02x'%mm)
    n2=r.i32()
    if r.ci and r.u8()!=1: raise ValueError('herm q')
    if n2!=n: raise ValueError('herm n2 %d'%n2)
    PD=r.fs(4*n); has2=r.i16(); D2=r.fs(2*n) if has2 else []
    t0,t1=r.f(),r.f(); tl=r.i16()
    if tl!=-2: raise ValueError('herm tail %d'%tl)
    p=I(deg,0,0,n,0)+D(*kn)+I(*mu)+I(n)+D(*PD)+h(has2)+D(*D2)+D(t0,t1)
  else:
    raise KeyError('class %d'%ot)
  if r.left(): raise ValueError('%d bytes left'%r.left())
  return ot,refs,p
