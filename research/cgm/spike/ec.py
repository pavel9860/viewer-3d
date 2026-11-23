import struct,bisect
import vi
D=lambda p,k,n=1: struct.unpack_from('<%dd'%n,p,k)
def parse_ec(refs,p):
  """edge curve (191/192/193): returns list of param maps: f_i(q) -> param on refs[i], q = param on refs[0]"""
  if (len(p)-3)%8 in (3,) and p[0]>=0x80: pass
  n,k,_=vi.rd(p,0)
  fl=struct.unpack_from('<H',p,k)[0]; k+=2
  a0,a1,tol=D(p,k,3); k+=24
  maps=[lambda q:q]
  for j in range(n):
    if fl&4:
      b0,b1=D(p,k,2); k+=16
      N,k,_=vi.rd(p,k); A=D(p,k,N); k+=8*N; B=D(p,k,N); k+=8*N
      def f(q,A=A,B=B):
        i=min(max(bisect.bisect_right(A,q)-1,0),len(A)-2)
        t=(q-A[i])/(A[i+1]-A[i]) if A[i+1]!=A[i] else 0
        return B[i]+t*(B[i+1]-B[i])
      maps.append(f)
    elif fl&2:
      b0,b1,s,off=D(p,k,4); k+=32
      maps.append(lambda q,s=s,off=off:(q-off)/s)
    elif fl&1:
      b0,b1,tol2=D(p,k,3); k+=24
      maps.append(lambda q:q)
    else: raise ValueError('ec flags %x'%fl)
  c0,c1,ctol=D(p,k,3); k+=24
  rest=p[k:]
  return maps,(a0,a1),rest
