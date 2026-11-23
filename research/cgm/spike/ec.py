import struct,bisect
import vi
D=lambda p,k,n=1: struct.unpack_from('<%dd'%n,p,k)
def count(p,k):
  # varint, or raw i32 in streams of level <= 16 (R8-R10, CXR6/7)
  if p[k]<0x80 and p[k+1:k+4]==b'\0\0\0' : return struct.unpack_from('<I',p,k)[0],k+4
  v,k,_=vi.rd(p,k); return v,k
def parse_ec(refs,p):
  """edge curve 191/192/193: maps f_i(q) -> param on refs[i] from q = param on refs[0]"""
  n,k=count(p,0)
  maps=[lambda q:q]; tail=False; ranges=[]
  for j in range(n):
    fl=struct.unpack_from('<H',p,k)[0]; k+=2
    a0,a1,tol=D(p,k,3); k+=24; ranges.append((a0,a1))
    if fl==0: continue
    if fl==1: b0,b1,_=D(p,k,3); k+=24; maps.append(lambda q:q); continue
    b0,b1=D(p,k,2); k+=16
    if fl==2: maps.append(lambda q,a0=a0,a1=a1,b0=b0,b1=b1: b0+(q-a0)*(b1-b0)/(a1-a0) if a1!=a0 else b0)
    elif fl==3:
      s,off=D(p,k,2); k+=16; maps.append(lambda q,s=s,off=off:(q-off)/s)
    elif fl==5:
      N,k=count(p,k); A=D(p,k,N); k+=8*N; B=D(p,k,N); k+=8*N
      def f(q,A=A,B=B):
        i=min(max(bisect.bisect_right(A,q)-1,0),len(A)-2)
        t=(q-A[i])/(A[i+1]-A[i]) if A[i+1]!=A[i] else 0
        return B[i]+t*(B[i+1]-B[i])
      maps.append(f)
    else: raise ValueError('ec flags %x'%fl)
    tail=True
  if tail: k+=24
  return maps,ranges,p[k:]
