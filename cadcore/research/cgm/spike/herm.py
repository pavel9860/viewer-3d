import struct
D=lambda p,k,n=1: struct.unpack_from('<%dd'%n,p,k)
I=lambda p,k: struct.unpack_from('<i',p,k)[0]
def parse_herm(p,k=0,dim=2):
  """piecewise Hermite: i32 deg, i32, i32, i32 nk, i32, knots[nk], mults i32[nk], i32 dim,
     P[dim][nk], D1[dim][nk], i16 has2, [D2[dim][nk]], range t0 t1"""
  deg=I(p,k); nk=I(p,k+12); kn=D(p,k+20,nk); k+=20+8*nk; m=[I(p,k+4*j) for j in range(nk)]; k+=4*nk
  cnt=I(p,k); k+=4
  assert cnt==nk, (cnt,nk)
  P=[D(p,k+8*nk*d,nk) for d in range(dim)]; k+=8*nk*dim
  D1=[D(p,k+8*nk*d,nk) for d in range(dim)]; k+=8*nk*dim
  has2=struct.unpack_from('<h',p,k)[0]; k+=2
  D2=None
  if has2: D2=[D(p,k+8*nk*d,nk) for d in range(dim)]; k+=8*nk*dim
  rng=D(p,k,2); k+=16
  return dict(deg=deg,knots=kn,mults=m,dim=dim,P=P,D1=D1,D2=D2,range=rng),k
def herm_eval(H,t):
  kn=H['knots']; n=len(kn)
  i=0
  while i<n-2 and t>kn[i+1]: i+=1
  h=kn[i+1]-kn[i]; s=(t-kn[i])/h if h else 0.0
  s2,s3,s4,s5=s*s,s**3,s**4,s**5
  if H['D2']:
    h0=1-10*s3+15*s4-6*s5; h1=s-6*s3+8*s4-3*s5; h2=0.5*s2-1.5*s3+1.5*s4-0.5*s5
    h5=10*s3-15*s4+6*s5; h4=-4*s3+7*s4-3*s5; h3=0.5*s3-s4+0.5*s5
    return [h0*H['P'][d][i]+h*h1*H['D1'][d][i]+h*h*h2*H['D2'][d][i]+h5*H['P'][d][i+1]+h*h4*H['D1'][d][i+1]+h*h*h3*H['D2'][d][i+1] for d in range(H['dim'])]
  a0=2*s3-3*s2+1; a1=s3-2*s2+s; b0=-2*s3+3*s2; b1=s3-s2
  return [a0*H['P'][d][i]+h*a1*H['D1'][d][i]+b0*H['P'][d][i+1]+h*b1*H['D1'][d][i+1] for d in range(H['dim'])]
