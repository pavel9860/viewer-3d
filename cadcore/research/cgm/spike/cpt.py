"""compact dialect (2011+) exploration: records with implicit ids, relative refs"""
import sys,os,glob,struct,collections
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__))); import new3,old
def sint(b,k):
  """signed compact int: low 2 bits 01 -> +(b>>2); 11 -> -((b>>2)+1); 00 -> (b>>2) bytes follow (LE
  magnitude, +); 10 -> same, negative"""
  x=b[k]; lo=x&3
  if lo==1: return x>>2,k+1
  if lo==3: return -((x>>2)+1),k+1
  n=x>>2; v=int.from_bytes(b[k+1:k+1+n],'little')
  if k+1+n>len(b): raise IndexError("cint past end")
  return (v if lo==0 else -v-1),k+1+n
def records(s):
  r,d,k=new3.chain(s)
  if r is None or k!=d[0]+d[1] or d[2]>82: return None
  out=[]; cur=0
  for (o,m,fl,t,ln,idf,b) in r:
    dv,_=sint(idf,0); cur+=dv
    out.append((cur,t,fl,s[b:b+ln]))
  return out,d
def files(level_max=82):
  for f in sorted(glob.glob(os.path.join(os.path.dirname(os.path.abspath(__file__)),'cgm2/*.cgm'))):
    s=open(f,'rb').read()
    if s.find(b'\xa8\x06\x7a')<0: continue
    try: x=records(s)
    except Exception: continue
    if x and x[1][2]<=level_max: yield f,s,x[0]
