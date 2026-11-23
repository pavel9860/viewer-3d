import struct
D=lambda p,k,n=1: struct.unpack_from('<%dd'%n,p,k)
I=lambda p,k: struct.unpack_from('<i',p,k)[0]
def kv(p,k):
  """knot vector: i32 deg, i32 a, i32 b, i32 nk, i32 c, f64 knots[nk], i32 mults[nk]"""
  deg=I(p,k); nk=I(p,k+12); kn=D(p,k+20,nk); k+=20+8*nk; m=[I(p,k+4*j) for j in range(nk)]; k+=4*nk
  full=[x for x,mm in zip(kn,m) for _ in range(mm)]
  return dict(deg=deg,knots=kn,mults=m,full=full,npoles=len(full)-deg-1),k
def surface(p):
  U,k=kv(p,0); V,k=kv(p,k)
  rat=struct.unpack_from('<h',p,k)[0]; k+=2
  nu,nv=U['npoles'],V['npoles']
  P=[D(p,k+24*(i*nv+j),3) for i in range(nu) for j in range(nv)]; k+=24*nu*nv
  W=None
  if rat: W=D(p,k,nu*nv); k+=8*nu*nv
  return dict(U=U,V=V,P=P,W=W,nu=nu,nv=nv),k
def curve(p,dim=3):
  U,k=kv(p,0)
  rat=struct.unpack_from('<h',p,k)[0]; k+=2
  n=U['npoles']; P=[D(p,k+8*dim*i,dim) for i in range(n)]; k+=8*dim*n
  W=None
  if rat: W=D(p,k,n); k+=8*n
  return dict(U=U,P=P,W=W),k
def basis(full,deg,t):
  n=len(full)-deg-1
  t=min(max(t,full[deg]),full[n])
  span=deg
  while span<n-1 and t>=full[span+1]: span+=1
  N=[0.0]*(deg+1); N[0]=1.0; left=[0.0]*(deg+1); right=[0.0]*(deg+1)
  for j in range(1,deg+1):
    left[j]=t-full[span+1-j]; right[j]=full[span+j]-t; saved=0.0
    for r in range(j):
      den=right[r+1]+left[j-r]; tmp=N[r]/den if den else 0.0
      N[r]=saved+right[r+1]*tmp; saved=left[j-r]*tmp
    N[j]=saved
  return span,N
def surf_eval(S,u,v):
  du,dv=S['U']['deg'],S['V']['deg']
  su,Nu=basis(S['U']['full'],du,u); sv,Nv=basis(S['V']['full'],dv,v)
  acc=[0.0,0.0,0.0,0.0]
  for a in range(du+1):
    for b in range(dv+1):
      i,j=su-du+a,sv-dv+b; w=S['W'][i*S['nv']+j] if S['W'] else 1.0; c=Nu[a]*Nv[b]*w; P=S['P'][i*S['nv']+j]
      acc[0]+=c*P[0]; acc[1]+=c*P[1]; acc[2]+=c*P[2]; acc[3]+=c
  return [acc[0]/acc[3],acc[1]/acc[3],acc[2]/acc[3]]
